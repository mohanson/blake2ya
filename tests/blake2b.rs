#[test]
fn blake2b_empty() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    let mut h = blake2::blake2b(p);
    let r = h.digest();
    let e = [
        0x78, 0x6a, 0x02, 0xf7, 0x42, 0x01, 0x59, 0x03, 0xc6, 0xc6, 0xfd, 0x85, 0x25, 0x52, 0xd2, 0x72, 0x91, 0x2f,
        0x47, 0x40, 0xe1, 0x58, 0x47, 0x61, 0x8a, 0x86, 0xe2, 0x17, 0xf7, 0x1f, 0x54, 0x19, 0xd2, 0x5e, 0x10, 0x31,
        0xaf, 0xee, 0x58, 0x53, 0x13, 0x89, 0x64, 0x44, 0x93, 0x4e, 0xb0, 0x4b, 0x90, 0x3a, 0x68, 0x5b, 0x14, 0x48,
        0xb7, 0x55, 0xd5, 0x6f, 0x70, 0x1a, 0xfe, 0x9b, 0xe2, 0xce,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_abc() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    let mut h = blake2::blake2b(p);
    h.update(b"abc");
    let r = h.digest();
    let e = [
        0xba, 0x80, 0xa5, 0x3f, 0x98, 0x1c, 0x4d, 0x0d, 0x6a, 0x27, 0x97, 0xb6, 0x9f, 0x12, 0xf6, 0xe9, 0x4c, 0x21,
        0x2f, 0x14, 0x68, 0x5a, 0xc4, 0xb7, 0x4b, 0x12, 0xbb, 0x6f, 0xdb, 0xff, 0xa2, 0xd1, 0x7d, 0x87, 0xc5, 0x39,
        0x2a, 0xab, 0x79, 0x2d, 0xc2, 0x52, 0xd5, 0xde, 0x45, 0x33, 0xcc, 0x95, 0x18, 0xd3, 0x8a, 0xa8, 0xdb, 0xf1,
        0x92, 0x5a, 0xb9, 0x23, 0x86, 0xed, 0xd4, 0x00, 0x99, 0x23,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_update() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    let mut h = blake2::blake2b(p);
    for _ in 0..128 {
        h.update(&[0; 64])
    }
    h.update(&[0; 32]);
    let r = h.digest();
    let e = [
        0x01, 0x7d, 0x62, 0x9e, 0x94, 0x59, 0xc2, 0x21, 0x32, 0x53, 0x35, 0xc8, 0x8c, 0xff, 0xca, 0xd3, 0xe5, 0x8c,
        0xbf, 0xaf, 0x08, 0x94, 0x36, 0xcb, 0x05, 0x04, 0x61, 0x1a, 0xe6, 0x54, 0x0b, 0x2b, 0x9e, 0x64, 0x72, 0x1a,
        0x21, 0xf0, 0x24, 0xe1, 0xac, 0x87, 0xf1, 0xdf, 0x63, 0x02, 0xb9, 0xe6, 0xc7, 0xa6, 0xbd, 0xa5, 0xb4, 0x7d,
        0x16, 0xd7, 0x04, 0xe9, 0x2a, 0x29, 0x32, 0x65, 0xe0, 0xfc,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_big() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    let mut h = blake2::blake2b(p);
    h.update(&[0; 8224]);
    let r = h.digest();
    let e = [
        0x01, 0x7d, 0x62, 0x9e, 0x94, 0x59, 0xc2, 0x21, 0x32, 0x53, 0x35, 0xc8, 0x8c, 0xff, 0xca, 0xd3, 0xe5, 0x8c,
        0xbf, 0xaf, 0x08, 0x94, 0x36, 0xcb, 0x05, 0x04, 0x61, 0x1a, 0xe6, 0x54, 0x0b, 0x2b, 0x9e, 0x64, 0x72, 0x1a,
        0x21, 0xf0, 0x24, 0xe1, 0xac, 0x87, 0xf1, 0xdf, 0x63, 0x02, 0xb9, 0xe6, 0xc7, 0xa6, 0xbd, 0xa5, 0xb4, 0x7d,
        0x16, 0xd7, 0x04, 0xe9, 0x2a, 0x29, 0x32, 0x65, 0xe0, 0xfc,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_digest_32() {
    let mut p = blake2::blake2b_params();
    p.digest(32);
    let mut h = blake2::blake2b(p);
    let r = &h.digest()[..32];
    let e = [
        0x0e, 0x57, 0x51, 0xc0, 0x26, 0xe5, 0x43, 0xb2, 0xe8, 0xab, 0x2e, 0xb0, 0x60, 0x99, 0xda, 0xa1, 0xd1, 0xe5,
        0xdf, 0x47, 0x77, 0x8f, 0x77, 0x87, 0xfa, 0xab, 0x45, 0xcd, 0xf1, 0x2f, 0xe3, 0xa8,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_key() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    p.key(&[
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11,
        0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x22, 0x23,
        0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
        0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
    ]);
    let mut h = blake2::blake2b(p);
    h.update(&[
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11,
        0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x22, 0x23,
        0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f, 0x30, 0x31, 0x32, 0x33, 0x34, 0x35,
        0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47,
        0x48, 0x49, 0x4a, 0x4b, 0x4c, 0x4d, 0x4e, 0x4f, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59,
        0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6a, 0x6b,
        0x6c, 0x6d, 0x6e, 0x6f, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7a, 0x7b, 0x7c, 0x7d,
        0x7e, 0x7f, 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
        0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f, 0xa0, 0xa1,
        0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf, 0xb0, 0xb1, 0xb2, 0xb3,
        0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf, 0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5,
        0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf, 0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7,
        0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf, 0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9,
        0xea, 0xeb, 0xec, 0xed, 0xee, 0xef, 0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb,
        0xfc, 0xfd, 0xfe,
    ]);
    let r = h.digest();
    let e = [
        0x14, 0x27, 0x09, 0xd6, 0x2e, 0x28, 0xfc, 0xcc, 0xd0, 0xaf, 0x97, 0xfa, 0xd0, 0xf8, 0x46, 0x5b, 0x97, 0x1e,
        0x82, 0x20, 0x1d, 0xc5, 0x10, 0x70, 0xfa, 0xa0, 0x37, 0x2a, 0xa4, 0x3e, 0x92, 0x48, 0x4b, 0xe1, 0xc1, 0xe7,
        0x3b, 0xa1, 0x09, 0x06, 0xd5, 0xd1, 0x85, 0x3d, 0xb6, 0xa4, 0x10, 0x6e, 0x0a, 0x7b, 0xf9, 0x80, 0x0d, 0x37,
        0x3d, 0x6d, 0xee, 0x2d, 0x46, 0xd6, 0x2e, 0xf2, 0xa4, 0x61,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_salt() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    p.salt(&[0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f]);
    let mut h = blake2::blake2b(p);
    h.update(b"abc");
    let r = h.digest();
    let e = [
        0x02, 0x6d, 0x34, 0x89, 0x6f, 0x69, 0x1f, 0xd4, 0xe5, 0x57, 0x76, 0x18, 0xf5, 0xa7, 0x11, 0x93, 0xcb, 0x3e,
        0xd1, 0xc9, 0xdf, 0x63, 0xba, 0x2c, 0x68, 0xcf, 0x65, 0x13, 0xf0, 0xd6, 0xe8, 0x31, 0x1d, 0x38, 0x32, 0xd9,
        0x4f, 0x4f, 0xd1, 0xad, 0xe2, 0x93, 0x6f, 0x08, 0x74, 0x05, 0xef, 0xaf, 0x91, 0x06, 0x9d, 0xdb, 0x89, 0x23,
        0x0f, 0x80, 0xa5, 0x95, 0x81, 0x06, 0xe7, 0x4c, 0x86, 0xc8,
    ];
    assert_eq!(r, e);
}

#[test]
fn blake2b_person() {
    let mut p = blake2::blake2b_params();
    p.digest(64);
    p.person(&[0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f]);
    let mut h = blake2::blake2b(p);
    h.update(b"abc");
    let r = h.digest();
    let e = [
        0x7b, 0xd5, 0xc1, 0xfc, 0xac, 0xfe, 0xd1, 0xe0, 0xf1, 0xe7, 0x30, 0x8a, 0xa4, 0x8c, 0xdb, 0xa7, 0x4b, 0x75,
        0xd1, 0xf6, 0xd4, 0x9d, 0x32, 0x37, 0x17, 0xb0, 0x16, 0x37, 0x29, 0xc1, 0x0d, 0x73, 0x03, 0x5b, 0x86, 0x59,
        0x2b, 0xe4, 0x4b, 0x07, 0x8d, 0x72, 0x95, 0x02, 0x49, 0x8e, 0x07, 0x3e, 0xb9, 0xe0, 0x71, 0x9a, 0x4b, 0x4a,
        0x12, 0x00, 0xda, 0x02, 0xb3, 0x7e, 0x8c, 0x34, 0xf3, 0x7a,
    ];
    assert_eq!(r, e);
}