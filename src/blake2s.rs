///  Message word schedule permutations for each round of both BLAKE2b and BLAKE2s are defined by SIGMA.
const BLAKE2S_SIGMA: [[u8; 16]; 10] = [
    [0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7, 0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf],
    [0xe, 0xa, 0x4, 0x8, 0x9, 0xf, 0xd, 0x6, 0x1, 0xc, 0x0, 0x2, 0xb, 0x7, 0x5, 0x3],
    [0xb, 0x8, 0xc, 0x0, 0x5, 0x2, 0xf, 0xd, 0xa, 0xe, 0x3, 0x6, 0x7, 0x1, 0x9, 0x4],
    [0x7, 0x9, 0x3, 0x1, 0xd, 0xc, 0xb, 0xe, 0x2, 0x6, 0x5, 0xa, 0x4, 0x0, 0xf, 0x8],
    [0x9, 0x0, 0x5, 0x7, 0x2, 0x4, 0xa, 0xf, 0xe, 0x1, 0xb, 0xc, 0x6, 0x8, 0x3, 0xd],
    [0x2, 0xc, 0x6, 0xa, 0x0, 0xb, 0x8, 0x3, 0x4, 0xd, 0x7, 0x5, 0xf, 0xe, 0x1, 0x9],
    [0xc, 0x5, 0x1, 0xf, 0xe, 0xd, 0x4, 0xa, 0x0, 0x7, 0x6, 0x3, 0x9, 0x2, 0x8, 0xb],
    [0xd, 0xb, 0x7, 0xe, 0xc, 0x1, 0x3, 0x9, 0x5, 0x0, 0xf, 0x4, 0x8, 0x6, 0x2, 0xa],
    [0x6, 0xf, 0xe, 0x9, 0xb, 0x3, 0x0, 0x8, 0xc, 0x2, 0xd, 0x7, 0x1, 0x4, 0xa, 0x5],
    [0xa, 0x2, 0x8, 0x4, 0x7, 0x6, 0x1, 0x5, 0xf, 0xb, 0x9, 0xe, 0x3, 0xc, 0xd, 0x0],
];

/// The initialization vector constant.
/// IV[i] = floor(2**w * frac(sqrt(prime(i+1)))), where prime(i) is the i:th prime number (2, 3, 5, 7, 11, 13, 17, 19).
const BLAKE2S_IV: [u32; 8] =
    [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];

/// Block bytes.
const BLAKE2S_BB: usize = 64;

/// G rotation constants.
const BLAKE2S_R1: u32 = 16;
const BLAKE2S_R2: u32 = 12;
const BLAKE2S_R3: u32 = 8;
const BLAKE2S_R4: u32 = 7;

/// The G primitive function mixes two input words, "x" and "y", into four words indexed by "a", "b", "c", and "d" in
/// the working vector v[0..15].
fn mixing(v: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize, x: u32, y: u32) {
    v[a] = v[a].wrapping_add(v[b]).wrapping_add(x);
    v[d] = (v[d] ^ v[a]).rotate_right(BLAKE2S_R1);
    v[c] = v[c].wrapping_add(v[d]);
    v[b] = (v[b] ^ v[c]).rotate_right(BLAKE2S_R2);
    v[a] = v[a].wrapping_add(v[b]).wrapping_add(y);
    v[d] = (v[d] ^ v[a]).rotate_right(BLAKE2S_R3);
    v[c] = v[c].wrapping_add(v[d]);
    v[b] = (v[b] ^ v[c]).rotate_right(BLAKE2S_R4);
}

/// Compression function F takes as an argument the state vector "h", message block vector "m" (last block is padded
/// with zeros to full block size, if required), 2w-bit offset counter "t", and final block indicator flag "f".  Local
/// vector v[0..15] is used in processing. F returns a new state vector. The number of rounds, "r", is 12 for BLAKE2b
/// and 10 for BLAKE2s. Rounds are numbered from 0 to r - 1.
fn reduce(h: &mut [u32; 8], m: &[u32; 16], t: &[u32; 2], f: &[u32; 2]) {
    let mut v = [0x00; 16];
    v[0x00..0x08].copy_from_slice(h);
    v[0x08..0x10].copy_from_slice(&BLAKE2S_IV);
    v[0x0c] ^= t[0];
    v[0x0d] ^= t[1];
    v[0x0e] ^= f[0];
    v[0x0f] ^= f[1];
    for s in BLAKE2S_SIGMA {
        mixing(&mut v, 0x0, 0x4, 0x8, 0xc, m[s[0x0] as usize], m[s[0x1] as usize]);
        mixing(&mut v, 0x1, 0x5, 0x9, 0xd, m[s[0x2] as usize], m[s[0x3] as usize]);
        mixing(&mut v, 0x2, 0x6, 0xa, 0xe, m[s[0x4] as usize], m[s[0x5] as usize]);
        mixing(&mut v, 0x3, 0x7, 0xb, 0xf, m[s[0x6] as usize], m[s[0x7] as usize]);
        mixing(&mut v, 0x0, 0x5, 0xa, 0xf, m[s[0x8] as usize], m[s[0x9] as usize]);
        mixing(&mut v, 0x1, 0x6, 0xb, 0xc, m[s[0xa] as usize], m[s[0xb] as usize]);
        mixing(&mut v, 0x2, 0x7, 0x8, 0xd, m[s[0xc] as usize], m[s[0xd] as usize]);
        mixing(&mut v, 0x3, 0x4, 0x9, 0xe, m[s[0xe] as usize], m[s[0xf] as usize]);
    }
    for i in 0..8 {
        h[i] = h[i] ^ v[i] ^ v[i + 8]
    }
}

/// Add n to message byte offset.
fn incoff(t: &mut [u32; 2], n: u32) {
    t[0] = t[0].wrapping_add(n);
    t[1] = t[1].wrapping_add((t[0] < n) as u32);
}

/// BLAKE2s parameter block structure.
pub struct Param2s {
    buf: [u8; 32],
    key: [u8; 32],
}

impl Param2s {
    /// Set digest byte length (1byte): an integer in [1, 64] for BLAKE2b, in [1, 32] for BLAKE2s.
    pub fn digest(&mut self, n: u8) {
        assert!(1 <= n && n <= 32);
        self.buf[0x00] = n;
    }

    /// Set key. Key length in [0, 64] for BLAKE2b, in [0, 32] for BLAKE2s.
    pub fn key(&mut self, n: &[u8]) {
        assert!(n.len() <= 32);
        self.buf[0x01] = n.len() as u8;
        self.key.copy_from_slice(n);
    }

    /// Salt is an arbitrary string of 16 bytes for BLAKE2b, and 8 bytes for BLAKE2s.
    pub fn salt(&mut self, n: &[u8; 8]) {
        self.buf[0x10..0x18].copy_from_slice(n);
    }

    /// Set personalization. An arbitrary string of 16 bytes for BLAKE2b, and 8 bytes for BLAKE2s.
    pub fn person(&mut self, n: &[u8; 8]) {
        self.buf[0x18..0x20].copy_from_slice(n);
    }
}

/// A context for computing the BLAKE2s checksum.
pub struct Blake2s {
    /// Internal state of the hash.
    h: [u32; 8],
    /// Message byte offset at the end of the current block.
    t: [u32; 2],
    /// Flag indicating the last block.
    f: [u32; 2],
    /// Buffer.
    b: [u8; BLAKE2S_BB],
    /// Buffer length.
    l: usize,
}

impl Blake2s {
    /// Update this hash object's state with the provided data.
    pub fn update(&mut self, data: &[u8]) {
        let mut data = data;
        if self.l != 0 {
            if data.len() <= BLAKE2S_BB - self.l {
                self.b[self.l..].copy_from_slice(data);
                self.l += data.len();
                return;
            }
            self.b[self.l..].copy_from_slice(&data[..BLAKE2S_BB - self.l]);
            incoff(&mut self.t, BLAKE2S_BB as u32);
            reduce(&mut self.h, unsafe { self.b.align_to::<u32>().1.try_into().unwrap() }, &self.t, &self.f);
            data = &data[BLAKE2S_BB - self.l..];
        }
        for _ in 0..data.len() / BLAKE2S_BB {
            self.b.copy_from_slice(&data[..BLAKE2S_BB]);
            incoff(&mut self.t, BLAKE2S_BB as u32);
            reduce(&mut self.h, unsafe { self.b.align_to::<u32>().1.try_into().unwrap() }, &self.t, &self.f);
            data = &data[BLAKE2S_BB..];
        }
        self.l = data.len();
        self.b[..self.l].copy_from_slice(data);
    }

    /// Return the digest value.
    pub fn digest(&mut self) -> [u8; 32] {
        self.b[self.l..].fill(0);
        self.f[0] = u32::MAX;
        incoff(&mut self.t, self.l as u32);
        reduce(&mut self.h, unsafe { self.b.align_to::<u32>().1.try_into().unwrap() }, &self.t, &self.f);
        unsafe { self.h.align_to::<u8>() }.1.try_into().unwrap()
    }
}

/// Create the parameter block of BLAKE2s. All general parameters are supported.
pub fn blake2s_params() -> Param2s {
    let mut r = Param2s { buf: [0; 32], key: [0; 32] };
    r.buf[0x02] = 0x01;
    r.buf[0x03] = 0x01;
    r
}

/// Core hasher state of BLAKE2s.
pub fn blake2s(param2s: Param2s) -> Blake2s {
    let mut r = Blake2s { h: [0; 8], t: [0; 2], f: [0; 2], b: [0; 64], l: 0 };
    for i in 0..8 {
        r.h[i] = BLAKE2S_IV[i] ^ u32::from_le_bytes(param2s.buf[i * 4..i * 4 + 4].try_into().unwrap())
    }
    if param2s.buf[1] != 0 {
        r.b[..32].copy_from_slice(&param2s.key);
        incoff(&mut r.t, BLAKE2S_BB as u32);
        reduce(&mut r.h, unsafe { r.b.align_to::<u32>().1.try_into().unwrap() }, &r.t, &r.f);
    }
    r
}
