"use strict";

// GLOBAL CONSTANTS...

const MAXVALUE = 9007199254740991; // 2**53 - 1
const SIEVEBUFBYTESZ = 512 * 1024; // 2**19

// UTILITY ARRAYS AND FUNCTIONS...

// faster than bit-twidling for masking of cull bits...
const bitMsk8 = Uint8Array.of(1, 2, 4, 8, 16, 32, 64, 128);

// faster than bit-twidling for masking of cull bits...
const nbitMsk8 = Uint8Array.of(0xFE,0xFD,0xFB,0xF7,0xEF,0xDF,0xBF,0x7F);

// faster than bit-twidling for "or'ing-out" upper bits for uint64 counts...
const cntMsk32 =
  Uint32Array.of( 0x00000001, 0x00000003, 0x00000007, 0x0000000F
                , 0x0000001F, 0x0000003F, 0x0000007F, 0x000000FF
                , 0x000001FF, 0x000003FF, 0x000007FF, 0x00000FFF
                , 0x00001FFF, 0x00003FFF, 0x00007FFF, 0x0000FFFF
                , 0x0001FFFF, 0x0003FFFF, 0x0007FFFF, 0x000FFFFF
                , 0x001FFFFF, 0x003FFFFF, 0x007FFFFF, 0x00FFFFFF
                , 0x01FFFFFF, 0x03FFFFFF, 0x07FFFFFF, 0x0FFFFFFF
                , 0x1FFFFFFF, 0x3FFFFFFF, 0x7FFFFFFF, 0xFFFFFFFF );

// modular 2/3/5 wheel gaps; two wheels for overflow...
const whlGaps = Uint8Array.of(6, 4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 4, 6, 2);

// modular 2/3/5 wheel residuals...
const whlRsds = Uint8Array.of(1, 7, 11, 13, 17, 19, 23, 29);

// index reverse lookup of above for one wheel of full circumference size...
const whlNdxs = Int8Array.of( -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3
                            , 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7 );

// Lookup table to find the upward rounded wheel position given the wheel
// modulo index; this covers the whole wheel circumference of 30 for two
// wheel spans to handle overflow; used for cull start address calculation...
const whlRndups = Int8Array.of( // impossible to index with zero!
  1, 1, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 13, 13, 17,
  17, 17, 17, 19, 19, 23, 23, 23, 23, 29, 29, 29, 29, 29, 29,
  31, 31, 37, 37, 37, 37, 37, 37, 41, 41, 41, 41, 43, 43, 47, 47,
  47, 47, 49, 49, 53, 53, 53, 53, 59, 59, 59, 59, 59, 59, 61, 61);

// parallel Lookup table to the above for new culling start wheel index;
// at the cost of a single extra lookup; it avoids several operations...
const whlRndupNdxs = Uint8Array.of(
  0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4,
  4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7,
  0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4,
  4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 0 );

function toOddNdx(x) { return (x - 1) >>> 1; } // for 32-bit math

// as long as `alpha` is at least 10.1 for `limit` of `MAXVALUE`,
// `alpha` is about 17.0 for `MAXVALUE` so this works for 32-bit math!...
function toWhlNdx(x) {
  const xd = Math.trunc(x / 30) | 0;
  const xm = (  x - (xd * 30)) | 0;
  return ((xd << 3) + whlNdxs[xm]) >>> 0; } // wheel index from value!

function ndxToPrm(x) { // 32-bit prime from value!
  return ((((x >>> 3) * 30) >>> 0) + whlRsds[x & 7]) >>> 0; }

// BLOCK COUNTING FUNCTIONS...

// because there is no hardware popcount built into JavaScript...
function popCount32(x) { // x must be a 32-bit integer!
  const x2 = (x - ((x >>> 1) & 0x55555555)) >>> 0;
  const x4 = ((x2 & 0x33333333) + ((x2 >>> 2) & 0x33333333)) >>> 0;
  const x8 = ((x4 & 0x0F0F0F0F) + ((x4 >>> 4) & 0x0F0F0F0F)) >>> 0;
  return (x8 * 0x01010101) >>> 24; }

// accumulates the phi counts of the 64-bit view of the sieve buffer;
// computation is from the low wheel (byte) index...
function syncBlkCntsFromFor(lowwhli, blkcnts, buf32, uprlmt) {
  const rellmt = (toWhlNdx(uprlmt)) - (lowwhli << 3);
  const cntlmt = (rellmt + 32) >>> 5; // rnd up to uint32 size!
  const ndxlmt = Math.min(buf32.length, cntlmt); // safety termination!
  let cnt = 0 >>> 0;
  for (let i = 0; i < ndxlmt; ++i) blkcnts[i + 1] = cnt += popCount32(buf32[i]);
  blkcnts[0] = 0; // just to be sure!
  return cnt; } // in case the total count is useful!

// accumulates the phi counts of the sixteen 32-bit view of the sieve buffer;
// computation is from the low wheel (byte) index...
function syncBlkCnts16FromFor(lowwhli, blkcnts, chnkcnts, uprlmt) {
  const rellmt = (toWhlNdx(uprlmt)) - (lowwhli << 3);
  const cntlmt = (rellmt + 512) >> 9; // round up to 16 * uint32 size!
  const ndxlmt = Math.min(chnkcnts.length, cntlmt); // safety termination!
  const lstblkcntsi = blkcnts.length - 1;
  blkcnts[lstblkcntsi] = 0; // set total count in case not set!
  let cnt = 0 >>> 0;
  for (let i = 0; i < ndxlmt; ++i) blkcnts[i + 1] = cnt += chnkcnts[i];
  blkcnts[0] = 0; // just to be sure!
  if (lowwhli == 0) blkcnts[lstblkcntsi] = cnt; // for zeroth page!
  return cnt; } // in case the total count is useful!

// use up-to-date block counts and the sieve buffer for phi count to point
// assume that `pnt` is within the given buffer starting at `lowwhli`...
function getPntCnt(lowwhli, blkcnts, buf32, pnt) {
  const relpntndx = (toWhlNdx(pnt)) - (lowwhli << 3);
  if (relpntndx < 0) return 0; // when `pnt` equals `lowwhli` * 30!
  const wrdi = relpntndx >>> 5; // 32-bit word index!
  const cnt = popCount32(buf32[wrdi] & cntMsk32[relpntndx & 31]);
  const rslt = blkcnts[wrdi] + cnt;
  return blkcnts[wrdi] + cnt; }

// use up-to-date block counts and the sieve buffer for phi count to point;
// this version counts within a block of sixteen 32-bit words and
// assumes that `pnt` is within the given buffer starting at `lowwhli`...
function getPntCnt16(lowwhli, blkcnts, buf32, pnt) {
  const relpntndx = (toWhlNdx(pnt)) - (lowwhli << 3);
  if (relpntndx < 0) return 0; // when `pnt` equals `lowwhli` * 30!
  const lstwrdi = relpntndx >>> 5;
  const blkwrdi = (lstwrdi & -16) >>> 0; // 16 32-bit words per block count!
  let cnt = popCount32(buf32[lstwrdi] & cntMsk32[relpntndx & 31]);
  for (let i = blkwrdi; i < lstwrdi; ++i) cnt += popCount32(buf32[i]);
  return blkcnts[blkwrdi >>> 4] + cnt; } // 16 32-bit words per block count!

// SIEVING FUNCTIONS...

// cull all representations of multiples of `bp` from `buf` up to `whlsz` bytes;
// `buf` is Uint8Array; `bp`**2 is < `buf` range (number of bytes times 30)
// `cntlmt` and `whlsz` are 32-bit integers; `cntlmt` is counting limit;
// `chnkcnts` is the count set bits per block of bytes/words in `buf`...
function cullBufOfTo(lowwhli, buf, bp, chnkcnts, cntlmt, whlsz) {
  const bpndx = toWhlNdx(bp); const bpwhli = bpndx >>> 3;
  let bpwpi = bpndx & 7; const bprsd = whlRsds[bpwpi];
  // do page start cull value calculation here...
  let cullval = (bp * bp) >>> 0; const lowval = (lowwhli * 30) >>> 0;
  if (cullval >= lowval) cullval -= lowval;
  else {
    cullval = (lowval - cullval) % (bp * 30);
    if (cullval != 0) {      
      const wpndx = bprsd + Math.trunc((cullval + bp - 1) / bp);
      const owp = whlRndups[wpndx];
      bpwpi = whlRndupNdxs[wpndx]; // adjust advanced culled state!
      cullval = bp * (owp - bprsd) - cullval;
    }
  }
  if (cntlmt > 0) { // this is complex counting cull!
    // round up to 8*uint64 with maximum size of whlsz!...
    const cntsz = Math.min(whlsz, (((toWhlNdx(cntlmt) + 512) & -512) >>> 3)
                                       - lowwhli);
    const lstchnkcntsndx = (chnkcnts.length - 1) >>> 0;; let totcnt = 0 >>> 0;
    for (let i = 0; i < 8; i += 1) { // by eight residual bit planes!
      const culli = toWhlNdx(cullval); const shft = culli & 7;
      let ci = culli >>> 3; const nmsk = nbitMsk8[shft];
      // from `chnkcnts` subtract the number of culls that are effective...
      for (; ci < cntsz; ci += bp) { // chunk size 512 bits!
        const v = buf[ci]; const cntbit = (v >>> shft) & 1;
        totcnt += cntbit; // total count adjustment!
        chnkcnts[ci >>> 6] -= cntbit; buf[ci] = v & nmsk; }
      for (; ci < whlsz; ci += bp) buf[ci] &= nmsk;
      cullval += bp * whlGaps[bpwpi]; bpwpi += 1; }
    return totcnt; } // return the total number of effective culls!
  else { // simple cull by eight residual bit planes!
    for (let i = 0; i < 8; i += 1) {
      const culli = toWhlNdx(cullval);
      const nmsk = nbitMsk8[culli & 7];
      for (let ci = culli >>> 3; ci < whlsz; ci += bp) buf[ci] &= nmsk;
      cullval += bp * whlGaps[bpwpi]; bpwpi += 1; }
    return 0; } // return value should never be used as not partial sieve!
}

// sync `chnkcnts` to current state of `buf`...
function syncChnkCntsFrom(chnkcnts, buf) {
  const buf32 = new Uint32Array(buf.buffer); // 32-bit view of `buf`!
  let totcnt = 0 >>> 0; const chnkcntsz = chnkcnts.length;
  for (let i = 0; i < chnkcntsz; ++i) {
    const blki = i << 4; let cnt = 0 >>> 0; // group of 16 32-bit words!
    for (let j = 0; j < 16; ++j) cnt += popCount32(buf32[blki + j]);
    chnkcnts[i] = cnt; totcnt += cnt; }
  return totcnt; // return total count in case it's useful!
}

// array of cull buffer culled of the pre cull primes;
// this is done by culling over small "wheels" that repeat,
// then copying this culled "wheel" ranges over the sieve buffer;
// `ptrnsz0` is the size of the pre-existing pattern;
// `buf` and `whlprms` are uint8 arrays and
// `chnkcnts` is the count set bits per 64 uint8 block of bytes/words...
function fillBufWithWhlPtrn(buf, ptrnsz0, whlprms, chnkcnts) {
  if (whlprms.length < 1) return 1;
  let ptrnsz = ptrnsz0 >>> 0; const bufsz = buf.length;
  for (const pcbp of whlprms) {
    // first copy last wheel size `pcbp` times; one size first time!...
    for (let rpt = 1; rpt < pcbp; ++rpt) {
      const strt = rpt * ptrnsz;
      if (strt >= bufsz) break;
      const sz = Math.min(ptrnsz, bufsz - strt);
      buf.copyWithin(strt, 0, sz); }
    // then cull new size buffer with `pcbp`; limit pcbp * circprdct...
    const bpi = toWhlNdx(pcbp); buf[bpi >>> 3] &= nbitMsk8[bpi & 7]; // bull bp!
    ptrnsz *= pcbp; cullBufOfTo(0, buf, pcbp, chnkcnts, 0, ptrnsz);
  }
  if (whlprms.length != 1) // for small ranges fill/init `chnkcnts` twice...
    return ptrnsz; // done if not final pre cull prime and not first of 7/11!
  // in case `ptrnsz` < `sievesz`, fill the rest of the sieve!...
  for (let strti = ptrnsz; strti < bufsz; strti += ptrnsz) {
      const sz = Math.min(ptrnsz, bufsz - strti);
      buf.copyWithin(strti, 0, sz); }
  // initialize the `chnkcnts` array for block counting...
  syncChnkCntsFrom(chnkcnts, buf);
  return ptrnsz;
}

// fill function to fill page-segment sieve buffer from wheel pattern buffer...
function fillBufFrom(lowwhli, buf, whlptrnbytesz, whlptrnbuf, chnkcnts) {
  const bufsz = buf.length >>> 0;
  for (let bufi = 0; bufi < bufsz; bufi += 16384) { // copy in 16K chunks
    const ptrni = ((lowwhli + bufi) % whlptrnbytesz) >>> 0;
    const cpysz = Math.min(16384, bufsz - bufi); // for safety!
    buf.set(whlptrnbuf.subarray(ptrni, ptrni + cpysz), bufi); }
  // sync `chnkcnts` to current state of `buf`...
  return syncChnkCntsFrom(chnkcnts, buf); // return total count!
}

// counting function is a generator that yields progress percentages..
function *countPrimesTo(limit) {

  // VALIDATE INPUT WITH DEFAULTS FOR SMALL LIMITS...

  if (limit < 49) {
    if (limit < 9) return limit < 2 ? 0 : limit < 3 ? 1 : (limit + 1) >>> 1;
    // simple Legendre algorithm up to where Gourdon implementation works...
    if (limit < 25)
      return ((limit - 1) >>> 1) - (Math.trunc(limit / 3 - 1) >>> 1) + 1;
    return ((limit + 1) >>> 1) - (Math.trunc(limit / 3 + 1) >>> 1)
              - (Math.trunc(limit / 5 + 1) >>> 1)
              + (Math.trunc(limit / 3 / 5 + 1) >>> 1) + 2;
  }

  // FUNCTION CONSTANTS...

  const sqrtlmt = Math.sqrt(limit) >>> 0;
  const cbrtlmt = Math.cbrt(limit) >>> 0;
  const loglmt = Math.log(limit);
  // `alpha` is about 5.0 for 1e12, 6.25 for 1e13,
  //     8.0 for 1e14, 10.25 for 1e15, and 13.5 for 1e16...
  const alpha = limit <= 10**12 ? 4.991429
                  : -90.7514286 + 9.53276388 * loglmt 
                      - 0.332764923 * loglmt**2 + 0.00409565096 * loglmt**3;
  const y = Math.min(sqrtlmt, (alpha * cbrtlmt) >>> 0);
  const ysz = toOddNdx(y) + 1;
  const sqrty = Math.sqrt(y) >>> 0;
  const sqrtysz = toOddNdx(sqrty) + 1;
  const sqrtlmtdivy = Math.sqrt(limit / y) >>> 0;
  const sqrtsqrtlmt = Math.sqrt(sqrtlmt) >>> 0;
  const lmtdivy = Math.trunc(limit / y);
  const lmtstar = Math.max(sqrtsqrtlmt, (limit / (y * y)) >>> 0);
  // round up to even number of eight 8-byte "blocks" for counting!...
  const maxsieve = toWhlNdx(lmtdivy) + 512;
  const sievebitsz = maxsieve - (maxsieve % 512);
  const fullsievesz = sievebitsz / 8; // wheel 2/3/5 index!
  const sievesz = Math.min(SIEVEBUFBYTESZ, fullsievesz); // wheel 2/3/5 index!

  // INITIALIZE SIEVING AND COUNTING ARRAYS NEEDED...

  const isprms32 = new Uint32Array(sievesz / 4);
  const isprms = new Uint8Array(isprms32.buffer);
  isprms[0] = 0xFF; // only need to initialize first wheel byte!
  const chnkcnts = new Uint16Array(sievesz / 64); // 16 uint32 per block!
  // last word of `blkcnts` is total count, special for partial sieving as
  // it may be zero when counting over less than the full sieve buffer;
  // `blkcnts` is used as 16 * uint32 counts for partial sieving and
  // individualuint32 counts for regular non-partial sieving, which is why
  // top element may be zero when counting over less than full buffer...
  const blkcnts = new Uint32Array(sievesz / 4 + 1); // uint32 count per block!

  // NECESSARY SIZE OF PRE CULL PRIMES BASED ON `SIEVESZ` AND `LIMIT`...

  // limit pre-cull primes for small counting range limits...
  const preCullOddPrms = function() {
    const maxPreCullOddPrms = [ 7, 11, 13, 17, 19 ];
    const sqrndx = maxPreCullOddPrms.findIndex(p => !(p**2 <= limit));
    const prms0 = sqrndx === -1 ? maxPreCullOddPrms
                    : maxPreCullOddPrms.slice(0, sqrndx);; let prod = 1;
    const prds = prms0.map(p => prod *= p);
    const prdndx = prds.findIndex(prd => !(prd <= sievesz));
    return prdndx === -1 ? prms0 : prms0.slice(0, prdndx);
  }();
  const k = preCullOddPrms.length + 3;

  // THE ACTUAL WORK STARTS HERE!!!...

  // CREATE FACTORS ARRAY USED WHEN PARTIAL SIEVING...

  const fctrs = new Uint16Array(ysz).fill(0xFFFF);
  for (const bp of [...[3, 5], ...preCullOddPrms]) {
      fctrs[toOddNdx(bp)] = 0; // cull `bp` itself phi style!
      for (let fi = toOddNdx(bp * bp); fi < ysz; fi += bp) fctrs[fi] = 0; }
  for (let bpi = 1; bpi < sqrtysz; ++bpi) {
      if (fctrs[bpi] >= 0xFFFF) {
        const bp = bpi + bpi + 1; // cull or xor all multiples of `bp`...
        for (let ci = toOddNdx(bp) + bp; ci < ysz; ci += bp)
          fctrs[ci] = fctrs[ci] >= 0xFFFF ? bp : fctrs[ci] ^ 1;
        const sqr = bp * bp; // cull to make square free!...
        for (let sqri = toOddNdx(sqr); sqri < ysz; sqri += sqr) fctrs[sqri] = 0;
      }
  }
  // apply moebius function for upper primes!...
  for (let bpi = sqrtysz; bpi < ysz; ++bpi) {
      if (fctrs[bpi] >= 0xFFFF) {
        const bp = bpi + bpi + 1;
        for (let ci = toOddNdx(bp) + bp; ci < ysz; ci += bp) fctrs[ci] ^= 1; } }
  fctrs[0] = 0xFFFE; // the "one" value is special as not a prime!
//  console.log("fctrs:", fctrs);

  // CREATE PRIME COUNTS ARRAY TO Y FOR "A" CALCULATION FOR `LIMIT//P//Q<=Y`...

  let cnt = k; const piys = new Int32Array(ysz);
  for (let bpi = 0; bpi < ysz; ++bpi)
    piys[bpi] = cnt += fctrs[bpi] >= 0xFFFF ? 1 : 0;
//  console.log("piys:", piys);

  // CREATE PRIMES ARRAY UP TO Y...

  const yprms = new Uint32Array(piys[piys.length - 1] - k);
  for (let bpi = 0; bpi < ysz; ++bpi)
    if (fctrs[bpi] >= 0xFFFF) yprms[piys[bpi] - k - 1] = bpi + bpi + 1;
//  console.log("yprms:", yprms);

  // CONSTANTS BASED ON THE `PIYS` ARRAY...

  const piy = piys[piys.length - 1];
  const picbrtlmt = piys[toOddNdx(cbrtlmt)];
  const pisqrtlmtdivy = piys[toOddNdx(sqrtlmtdivy)];
  const pilmtstar = piys[toOddNdx(lmtstar)];
  const pisqrty = piys[toOddNdx(sqrty)];
  const picbrtlmtdivy = Math.max(k, piys[toOddNdx(Math.cbrt(limit / y))]);
//  console.log(piy, picbrtlmt, pisqrtlmtdivy, pilmtstar, pisqrty, picbrtlmtdivy);

  // PRE-SIEVE THE SIEVE BUFFER FOR PHI0 AND SYNCH WITH BLOCK COUNT ARRAY...

  const psz = fillBufWithWhlPtrn( isprms, 1
                                , preCullOddPrms.slice(0, -1), chnkcnts );
  syncBlkCntsFromFor(0, blkcnts, isprms32, psz * 30);
//  console.log(psz, isprms, isprms32, blkcnts);

  // CALCULATE GOURDON'S "PHI0" TERM USING NOW INITIALIZED VALUES/ARRAYS...

  const phi0 = function() {
    const lstpreprm = preCullOddPrms[preCullOddPrms.length - 1];
    const kprmscirc =
      preCullOddPrms.slice(0, -1).reduce((acc, p) => acc * p, 1) * 15;
    const kprmshits =
      preCullOddPrms.slice(0, -1).reduce((acc, p) => acc * (p - 1), 1) * 8;
    // tinyphi function uses constants from k determinations above...
    function tinyphi(x) {
      function tp(v) { // this could be well beyond 32-bit range!
        const oxndv = Math.trunc((v - 1) / 2);
        const circdiv = Math.trunc(oxndv / kprmscirc);
        const circmod = oxndv - circdiv * kprmscirc;
        return circdiv * kprmshits
                 + getPntCnt(0, blkcnts, isprms32, (circmod << 1) + 1); }
      return tp(x) - tp(Math.trunc(x / lstpreprm)); }
    // `phi0` can be calculated just by doing the moebius summation over the
    // limit divided by the numbers represented by qualified`fctrs`...
    let rslt = 0;
    for (let fctri = 0; fctri < ysz; ++fctri) { // scan all factors...
      const fctr = fctrs[fctri];
      if (fctr <= 1) continue; // "pk" factor or not square free: ignore!
      const mbv = ((fctr & 1) << 1) - 1; // inverse moebius value!
      rslt -= mbv * tinyphi(Math.trunc(limit / ((fctri << 1) + 1)));
    }
    return rslt;
  }();
//  console.log(phi0);

  // THE TIME CONSUMING WORK IS TO CALCULATE A, B, C, AND D TERMS...
  
  let aacc = 0; let bacc = 0; let cacc = 0; let dacc = 0;

  // INITIALIZE SIEVING/COUNTING/`WHLPTRN` ARRAYS AS NEEDED...

  // the wheel pattern is in "phi" form, so that it maps to every page...
  const whlptrnbytesz = psz * preCullOddPrms.slice(-1);
  const whlptrn = function() {    
    if (fullsievesz > sievesz) { // need preparation for page segmentation...
      const whlptrnsz = (whlptrnbytesz + 63 + 16384) & -64; // rnd 16 * uint32!
      const arr = new Uint8Array(whlptrnsz); // extra 16k for copy overflows!
      arr.set(isprms.subarray(0, psz), 0); // copy in existing pattern!
      // finish full wheel pattern and fill `isprms` array with it...
      fillBufWithWhlPtrn(arr, psz, preCullOddPrms.slice(-1), chnkcnts);
      fillBufFrom(0, isprms, whlptrnbytesz, arr, chnkcnts);
//    console.log(getPntCnt(0, blkcnts, isprms32, totptrnsz * 30));
      return arr;
    }
    else {
      fillBufWithWhlPtrn(isprms, psz, preCullOddPrms.slice(-1), chnkcnts);
//    console.log(isprms32, chnkcnts, blkcnts);
      return isprms.subarray(); // used by `fillBufWithWhlPtrn` for `revBPntGen`
    }
  }();
//  console.log(getPntCnt16(0, blkcnts, isprms32, (totpsz + 16384) * 30));

  // CREATE A REVERSE POINT GENERATOR FOR SQRT(`LIMIT`) >= PRIME > `Y` WHERE
  // POINT IS `LIMIT` DIVIDED BY PRIME FOR THE "B" CALCULATION...

  // works without "pi" form as `Y` is always higher than highest "pk"...
  let numrevbpnts = 0; // accumulates number of reverse points over range
  const revBPntGen = function*() {
    const hindx = toWhlNdx(sqrtlmt);
    const bufhisz = ((hindx + 512) & -512) >>> 3; // round up to nearest block!
    const buflowhli = sievesz >= fullsievesz // multi-page: round down to block!
                        ? 0 : (toWhlNdx(y) & -512) >>> 3;
    const svrng = bufhisz - buflowhli;
    const bufsz = Math.min(SIEVEBUFBYTESZ, svrng);
    const numbufpgs = Math.ceil(svrng / bufsz); // rounded up!
    const strtwhli = buflowhli + (numbufpgs - 1) * bufsz; // from top page!
    const histrtbpi = hindx - (strtwhli << 3);
    const buf32 = new Uint32Array(bufsz / 4);
    const buf = new Uint8Array(buf32.buffer);
    const bufchnkcnts = new Uint16Array(bufsz >>> 6); // 16 uint32's per block
    for (let lwwhli = strtwhli; lwwhli >= buflowhli; lwwhli -= bufsz) {
      const frstval = lwwhli * 30;
      const nxtval = lwwhli >= strtwhli ? bufhisz * 30 : frstval + bufsz * 30;
      const bpisz = piys[toOddNdx(Math.sqrt(nxtval - 1))] - k;
      const totcnt = fillBufFrom( lwwhli, buf, whlptrnbytesz
                                , whlptrn, bufchnkcnts );
      const cullsz = lwwhli >= strtwhli ? bufhisz - lwwhli : bufsz;
      for (let bpi = 0; bpi < bpisz; ++bpi) {
        const bp = yprms[bpi];
        cullBufOfTo(lwwhli, buf, bp, bufchnkcnts, 0, cullsz); }
      // loop over indices from rep for sqrtlmt down to y yielding "B" points...
      const strtbpi = lwwhli >= strtwhli ? histrtbpi : (bufsz << 3) - 1;
      const bsprm = lwwhli * 30;
      // this generator iteration is a bit slow but acceptable...
      for (let bpi = strtbpi; bpi >= 0; --bpi) {
        if ((buf[bpi >>> 3] >>> (bpi & 7)) & 1) {
          const bp = bsprm + ndxToPrm(bpi); if (bp <= y) return;
          numrevbpnts += 1; yield Math.trunc(limit / bp); }
      }
    }
    return;
  }(); // create generator and...
  let revBPnt = revBPntGen.next(); // advance to first value or "done"!

  // COMPUTE BY PAGE SEGMENTS OVER THE RANGE OF PAGES STARTING AT `CBRTLMT`...

  // initialize base pi's array for each `yprm` up to `lmtstar` used for
  // "D" accumulation calculations, with extra element for total phi;
  // these are updated per segment page calculation...
  const lstbasephisi = pilmtstar - k;
  const basephis = new Int32Array(lstbasephisi + 1);
  basephis[lstbasephisi] = k - 1; // last element set to "pi"!

/*
  // `numpgs`, `pgnum`, and `pgnumchk` are used for progress computation...
  const numpgs = Math.ceil(fullsievesz / SIEVEBUFBYTESZ);
  let pgnum = 0; let pgnumchk = 10; // update progress every 10 pages!
// */
  // `aprxnumops`, `numops`, and `numopschk` are used for progress computation...
  const aprxnumops = function() { // slightly overestimate number of ops...
    const numprmstocbrt = cbrtlmt / Math.log(cbrtlmt);
    return Math.trunc((numprmstocbrt * (numprmstocbrt - 1)) / 2 * 4.2);
  }();
  let numops = 0; let numopschk = 0; // update progress every 10**7 ops!
  for (let lwwhli = 0; lwwhli < fullsievesz; lwwhli += sievesz) {

    // CHECK PROGRESS ONCE PER PAGE PLUS EVERY BASE PRIME FOR "A" CALCULATION
    // WHERE COUNTS ARE MUCH DENSER THAN ANYWHERE ELSE...

    if (numops >= numopschk) {
      numopschk += 100000000; yield numops / aprxnumops * 100; }

    // CONSTANTS NEEDED FOR "D1", "D2","C1", "C2", AND "A" CALCULATIONS...
    
    const lstblkcnti = blkcnts.length - 1;
    const frstval = lwwhli * 30; const nxtval = frstval + sievesz * 30;
    const grddfrstval = Math.max(frstval, 1);
    const grddnxtval = Math.min(nxtval, lmtdivy);
    const maxprmsqrd = Math.trunc(limit / grddfrstval);
    const maxprm = Math.trunc(Math.sqrt(maxprmsqrd));
    const acsplit = Math.min(nxtval, sqrtlmt);
    const lmtdivacsplit = Math.trunc(limit / acsplit);
    const pisqrtfrstval = piys[toOddNdx(Math.trunc(Math.sqrt(grddfrstval)))];
    const pimaxc2 = piys[toOddNdx(Math.min(lmtdivacsplit / y), lmtstar)];

    // PAGE SEGMENT PARTIAL CULLING AND "D" COMPUTATION STARTS HERE...

    // page segment already filled with pattern for zeroth page;
    // fill with wheel pattern for all subsequent pages!...
    let curphi = function() {
      if (lwwhli <= 0) {
        syncBlkCnts16FromFor(0, blkcnts, chnkcnts, sievesz * 30);
        return blkcnts[lstblkcnti]; }
      const picnt =
        fillBufFrom(lwwhli, isprms, whlptrnbytesz, whlptrn, chnkcnts);
      syncBlkCnts16FromFor(lwwhli, blkcnts, chnkcnts, nxtval);
      return picnt;
    }();
    const bpisz = piys[toOddNdx(Math.sqrt(Math.min(nxtval - 1, lmtdivy)))] - k;
    for (let bpi = 0; bpi < bpisz; ++bpi) {
      const bp = yprms[bpi];

      // constants used for both "D1" and "D2" calculations...
      const lmtdivbplow = Math.min(Math.trunc(maxprmsqrd / bp), y);
      const lmtdivbphigh = Math.min(Math.trunc(limit / (bp * grddnxtval)), y);
      const maxm = Math.min(Math.trunc(limit / bp**3), lmtdivbplow);

      if (bp <= sqrty) { // COMPUTE "D1" FOR BASE PRIMES <= ISQRT(Y)...
        const maxfctri = toOddNdx(maxm);
        const minfctri = toOddNdx(Math.max(lmtdivbphigh, Math.trunc(y / bp)));
        const totprvphi = basephis[bpi] - bpi; // remove count for base primes!
//*
        for (let fctri = maxfctri; fctri > minfctri; --fctri) {
          const fctr = fctrs[fctri];
          if (fctr > bp) {
            const mbv = ((fctr & 1) << 1) - 1; // moebius value for this factor!!
            const bpm = bp * ((fctri << 1) + 1);
            dacc += mbv * (getPntCnt16( lwwhli, blkcnts, isprms32
                                      , Math.trunc(limit / bpm) ) + totprvphi);
            numops += 1; } }
// */
        basephis[bpi] += curphi; // add phi for partial sieve for base prime!
        curphi -= cullBufOfTo( lwwhli, isprms, bp
                             , chnkcnts, nxtval, sievesz );
        syncBlkCnts16FromFor(lwwhli, blkcnts, chnkcnts, nxtval);
      } else {
        if (bp <= lmtstar) {
          // ADD "D2" CONTRIBUTION FOR ABOVE < BASE PRIMES <= `LIMIT`^(1/4);
          // IS PHI'S DUE TO UNIQUE PRIME PAIRS (UPP) - FROM YPRMS ARRAY...
          const maxp2i = piys[toOddNdx(maxm)] - k - 1;
          const minp2i = piys[toOddNdx(Math.min( Math.max(lmtdivbphigh, bp)
                                               , maxm ))] - k;
          const totprvphi = basephis[bpi] - bpi; // remove cnt for base primes!
          numops += maxp2i - minp2i + 1;
          for (let p2i = maxp2i; p2i >= minp2i; --p2i) {
            const bpm = bp * yprms[p2i];
//*
            dacc += totprvphi + getPntCnt16( lwwhli, blkcnts, isprms32
                                           , Math.trunc(limit / bpm) );
// */
          }
          basephis[bpi] += curphi; // add phi for partial sieve for base prime!
          curphi -= cullBufOfTo(lwwhli, isprms, bp, chnkcnts, nxtval, sievesz);
          const upprlmt = Math.trunc(limit / bp ** 2);
          syncBlkCnts16FromFor(lwwhli, blkcnts, chnkcnts, upprlmt);
        } else // regular sieve for the rest of the bp's...
          cullBufOfTo(lwwhli, isprms, bp, chnkcnts, nxtval, sievesz);
      }
    }

    // ADJUST BASEPHIS FOR ELEMENTS BETWEEN INDEX BPISZ AND PILMTSTAR...

    for (let bpi = bpisz; bpi < pilmtstar - k; ++bpi) basephis[bpi] += curphi;

    // FINAL SYNCH OF BLOCK COUNTS AND SET TOTAL PREVIOUS PI...

    syncBlkCntsFromFor(lwwhli, blkcnts, isprms32, nxtval);
    // pick up the total pi count from the total phi count;
    // upper element of `basephis` is total previous pi not phi as set below...
    const totprvpi = basephis[lstbasephisi];
    // SIEVE BUFFER NOW FULLY CULLED/SIEVED AND BLOCK COUNTS SYNC'ED FOR "PI"!

    // COMPUTE "C1" CONTRIBUTION FOR BASE PRIMES <= ISQRT(Y) USING ABOVE ARRAYS;
    // "C1", "C2", AND "A" ARE PARTS TO ONLY SIEVE TO LIMIT / `LMTSTAR`**2...

    for (let bpi = picbrtlmtdivy - k; bpi < pisqrty - k; ++bpi) {
      const bp = yprms[bpi];
      const maxm = Math.trunc(Math.min(maxprmsqrd / bp, limit / bp**2, y));
      const minm = Math.trunc(Math.max(limit / bp**3, y / bp));
      const maxfctri = toOddNdx(maxm); const minfctri = toOddNdx(minm);
      for (let fctri = maxfctri; fctri > minfctri; --fctri) {
        const fctr = fctrs[fctri];
        if (fctr > bp) {
          const mbv = ((fctr & 1) << 1) - 1; // moebius value for this factor!!
          const bpm = bp * ((fctri << 1) + 1);
//*
          cacc += mbv * ( getPntCnt( lwwhli, blkcnts, isprms32
                                   , Math.trunc(limit / bpm) )
                            + totprvpi - bpi - k + 1 ); // xtra 1 already added!
// */
          numops += 1; } }
    }

    // ADD "C2" CONTRIBUTION FOR ABOVE < BASE PRIMES <= ISQRT(ISQRT(LIMIT)); WHICH
    // IS PI CONTRIBUTIONS DUE UNIQUE PRIME PAIRS (UPP) - FROM `YPRMS` ARRAY USING
    // FULLY CULLED SIEVE BUFFER AND FINALIZED BLOCK COUNTING ARRAYS...

    const maxp1isz = piys[toOddNdx(Math.min(maxprm, lmtstar))] - k;
    const minp1i = Math.max(picbrtlmtdivy, pisqrty, pisqrtfrstval, pimaxc2)- k;
    for (let p1i = minp1i; p1i < maxp1isz; ++p1i) {
      if ((p1i & 511) == 0 && numops >= numopschk) { // progress indication...
        numopschk += 1000000; yield numops / aprxnumops * 100; }

      const p1 = yprms[p1i];
      const maxp2 = Math.trunc(Math.min(maxprmsqrd / p1, limit / p1**2, y));
      const maxp2i = piys[toOddNdx(maxp2)] - k - 1;
      const minp2 = Math.trunc(Math.min(Math.max( lmtdivacsplit / p1
                                                , limit / p1**3, p1 ), maxp2));
      const minp2i = piys[Math.min(ysz - 1, toOddNdx(minp2))] - k;
      numops += maxp2i - minp2i + 1;
      for (let p2i = maxp2i; p2i >= minp2i; --p2i) {
        const dvsr = p1 * yprms[p2i];
//*
        cacc += getPntCnt(lwwhli, blkcnts, isprms32, Math.trunc(limit / dvsr))
                  + totprvpi - p1i - k + 1; // extra count already added!
// */
      }
    }

    // CALCULATE "A" USING FULLY CULLED SIEVE BUFFER AND SYNCH'ED BLOCK COUNTS...
    // "A" is completely calculated when pages reach `sqrtlmt`...
    if (frstval <= sqrtlmt) {
      const minp1i =
        piys[ toOddNdx( Math.max( lmtstar
                                , Math.min( Math.trunc(lmtdivacsplit / acsplit)
                                          , cbrtlmt ) ) ) ] - k;
      const p1isz = piys[toOddNdx(Math.min(cbrtlmt, maxprm))] - k;
      for (let p1i = minp1i; p1i < p1isz; ++p1i) {
        if ((p1i & 15) == 0 && numops >= numopschk) { // progress indication...
          numopschk += 20000000; yield numops / aprxnumops * 100; }

          const p1 = yprms[p1i]; const lmtdivp1 = Math.trunc(limit / p1);
        const sqrtlmtdivp1 = Math.trunc(Math.sqrt(lmtdivp1));
        const maxp2 = Math.min(sqrtlmtdivp1, Math.trunc(maxprmsqrd / p1));
        const minp2 = Math.min(Math.trunc(lmtdivacsplit / p1), sqrtlmtdivp1);
        const maxp2i = piys[toOddNdx(maxp2)] - k - 1;
        const minp2i = piys[toOddNdx(Math.max(minp2, p1))] - k;
        const splitsz = piys[toOddNdx(Math.max(lmtdivp1 / y, minp2, p1))] - k;
        const aftspliti = Math.min(splitsz - 1, maxp2i);
//*
        // reverse order so points are increasing; better cache associativity...
        numops += maxp2i - minp2i + 1;
        for (let p2i = maxp2i; p2i >= splitsz; --p2i)
          aacc += 2*( totprvpi + getPntCnt( lwwhli, blkcnts, isprms32
                                          , Math.trunc(lmtdivp1 / yprms[p2i]) ) );
        for (let p2i = aftspliti; p2i >= minp2i; --p2i)
          aacc += totprvpi + getPntCnt( lwwhli, blkcnts, isprms32
                                      , Math.trunc(lmtdivp1 / yprms[p2i]) );
// */
      }
    }

    // CALCULATE "B" USING FULLY CULLED SIEVE BUFFER AND SYNCH'ED BLOCK COUNTS...

    let lstbpnt = revBPnt.value;
    for ( ; !revBPnt.done && revBPnt.value < nxtval;
          revBPnt = revBPntGen.next()) {
//*
      bacc += totprvpi + getPntCnt(lwwhli, blkcnts, isprms32, revBPnt.value);
// */
      lstbpnt = revBPnt.value; }

    // UPDATE THE TOTAL PHI IN THE `BASEPHIS` ARRAY (LAST)...

    basephis[lstbasephisi] += blkcnts[lstblkcnti];

  }
//  console.log(numpgs, pgnum);

  // CALCULATE GOURDON'S "SIGMA" TERM USING NOW INITIALIZED VALUES/ARRAYS...

  numops += numrevbpnts;
  const pisqrtlmt = piy + numrevbpnts;
//  console.log("new pisqrtlmt:", pisqrtlmt);
  const sigma = function(a, b, c, d, e) {
    const sigarr = new Float64Array(7);
    sigarr[0] = a - 1 + Math.trunc((e * (e - 1)) / 2)
                      - Math.trunc((a * (a - 1)) / 2);
    sigarr[1] = Math.trunc(((a - b) * (a - b - 1)) / 2);
    sigarr[2] = a * (b - c - Math.trunc((c * (c - 3)) / 2)
                           + Math.trunc((d * (d - 3)) / 2));
    sigarr[3] = Math.trunc((b * (b - 1) * (b + b - 1)) / 6) - b
                  - Math.trunc((d * (d - 1) * (d + d - 1)) / 6) + d;
    const strti = d - k; const spltsz = c - k; const endsz = b - k;
    for (let yi = strti; yi < spltsz; ++yi)
      sigarr[4] += piys[toOddNdx(Math.trunc(limit / (yprms[yi] * y)))];
    for (let yi = spltsz; yi < endsz; ++yi)
      sigarr[5] += piys[toOddNdx(Math.trunc(limit / (yprms[yi] ** 2)))];
    for (let yi = strti; yi < endsz; ++yi)
      sigarr[6] -=
        piys[toOddNdx(Math.trunc(Math.sqrt(limit / (yprms[yi]))))] ** 2;
    sigarr[4] *= a; return sigarr.reduce((acc, val) => acc + val, 0);
  }(piy, picbrtlmt, pisqrtlmtdivy, pilmtstar, pisqrtlmt);

  // ASSEMBLE FINAL RESULT...

  console.log(limit, sqrtlmt, cbrtlmt, alpha, y, ysz, sqrty, sqrtysz, sqrtsqrtlmt, sqrtlmtdivy, lmtdivy, lmtstar, fullsievesz, preCullOddPrms, k);
  console.log(piy, picbrtlmt, pisqrtlmtdivy, pilmtstar, pisqrtlmt);
  console.log(aacc, bacc, cacc, dacc, phi0, sigma);
  console.log("Total number of divisions:", numops);
  return aacc - bacc + cacc + dacc + phi0 + sigma; // classic Gourdon formula!
}

let cancelled = false;

function doit() {
  const limit =  Math.floor(parseFloat(document.getElementById('limit').value));
  const start = Date.now();
  const prgrsgen = countPrimesTo(limit);
  function prgrsfnc() {
    if (cancelled) {
      document.getElementById('output').innerText = "Cancelled!!!";
    }
    else {
      const prgrs = prgrsgen.next();
      if (!prgrs.done) { // pgrslt.value is a percent done...            
        document.getElementById('output').innerText =
          "Counted " + (prgrs.value.toFixed(3)) + "%";
        setTimeout(prgrsfnc, 0); return;
      }
      // done; value is the count result...
      const elpsd = Date.now() - start;
      document.getElementById('output').innerText =
        "Found " + prgrs.value + " primes up to "
          + limit + " in " + elpsd + " milliseconds.";
    }
    cancelled = false;    
    document.getElementById('go').onclick = strtclick;
    document.getElementById('go').value = "Start Counting...";            
    document.getElementById('go').disabled = false;
    return;
  }
  prgrsfnc();
}

const cancelclick = function () {
  cancelled = true;
  document.getElementById('go').disabled = true;
  document.getElementById('go').value = "Cancelled!!!";
  document.getElementById('go').onclick = strtclick;
}

const strtclick = function () {
  const limit =  Math.floor(parseFloat(document.getElementById('limit').value));
  if (!Number.isInteger(limit) || (limit < 0) || (limit > MAXVALUE)) {    
    document.getElementById('output').innerText =
      "Top limit must be an integer between 0 and " + MAXVALUE + "!";
    return;
  }
  document.getElementById('output').innerText = "Counted 0.000%";
  document.getElementById('go').onclick = cancelclick;
  document.getElementById('go').value = "Running, click to cancel...";
  cancelled = false;
  setTimeout(doit, 0);
};

document.getElementById('go').onclick = strtclick;
document.addEventListener("keydown", function(event) {
  // Check if Enter is pressed
  if (event.key === "Enter") {
    // Make sure the user isn't typing in a multi-line paragraph block (textarea)
    if (document.activeElement.tagName !== "TEXTAREA") {
      event.preventDefault();
      document.getElementById("go").click();
    }
  }
});
