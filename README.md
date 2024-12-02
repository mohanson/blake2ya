# Blake2 Yet Another

The [BLAKE2](https://www.rfc-editor.org/rfc/rfc7693.html) cryptographic hash and message authentication code.

```toml
[dependencies]
blake2ya = "0.9"
```

## Example

```rs
let mut p = blake2ya::blake2b_params();
p.digest(64);
let mut h = blake2ya::blake2b(p);
h.update(b"abc");
let mut r = [0; 64];
h.digest(&mut r);
let e = [
    0xba, 0x80, 0xa5, 0x3f, 0x98, 0x1c, 0x4d, 0x0d, 0x6a, 0x27, 0x97, 0xb6, 0x9f, 0x12, 0xf6, 0xe9, 0x4c, 0x21,
    0x2f, 0x14, 0x68, 0x5a, 0xc4, 0xb7, 0x4b, 0x12, 0xbb, 0x6f, 0xdb, 0xff, 0xa2, 0xd1, 0x7d, 0x87, 0xc5, 0x39,
    0x2a, 0xab, 0x79, 0x2d, 0xc2, 0x52, 0xd5, 0xde, 0x45, 0x33, 0xcc, 0x95, 0x18, 0xd3, 0x8a, 0xa8, 0xdb, 0xf1,
    0x92, 0x5a, 0xb9, 0x23, 0x86, 0xed, 0xd4, 0x00, 0x99, 0x23,
];
assert_eq!(r, e);
```
