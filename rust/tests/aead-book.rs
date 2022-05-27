use hacl_rust::aead::{self, Aead, Algorithm};
#[test]
fn stateful() {
    // ANCHOR: stateful
    let key = [
        0x5b, 0x96, 0x04, 0xfe, 0x14, 0xea, 0xdb, 0xa9, 0x31, 0xb0, 0xcc, 0xf3, 0x48, 0x43, 0xda,
        0xb9, 0x5b, 0x96, 0x04, 0xfe, 0x14, 0xea, 0xdb, 0xa9, 0x31, 0xb0, 0xcc, 0xf3, 0x48, 0x43,
        0xda, 0xb9,
    ];
    // ANCHOR: stateful_cipher
    let cipher = Aead::new(Algorithm::Chacha20Poly1305, &key).unwrap();
    // ANCHOR_END: stateful_cipher

    let iv = [
        0x02, 0x83, 0x18, 0xab, 0xc1, 0x82, 0x40, 0x29, 0x13, 0x81, 0x41, 0xa2,
    ];
    let msg = [
        0x00, 0x1d, 0x0c, 0x23, 0x12, 0x87, 0xc1, 0x18, 0x27, 0x84, 0x55, 0x4c, 0xa3, 0xa2, 0x19,
        0x08,
    ];
    let aad = [];

    // ANCHOR: stateful_encrypt
    let (ciphertext, tag) = cipher.encrypt(&msg, &iv, &aad).unwrap();
    let msg_ = cipher.decrypt(&ciphertext, &tag, &iv, &aad).unwrap();
    // ANCHOR_END: stateful_encrypt

    assert_eq!(&msg[..], &msg_[..]);
    // ANCHOR_END: stateful
}

#[test]
fn single_shot() {
    // ANCHOR: single_shot
    let key = [
        0x5b, 0x96, 0x04, 0xfe, 0x14, 0xea, 0xdb, 0xa9, 0x31, 0xb0, 0xcc, 0xf3, 0x48, 0x43, 0xda,
        0xb9, 0x5b, 0x96, 0x04, 0xfe, 0x14, 0xea, 0xdb, 0xa9, 0x31, 0xb0, 0xcc, 0xf3, 0x48, 0x43,
        0xda, 0xb9,
    ];
    let iv = [
        0x02, 0x83, 0x18, 0xab, 0xc1, 0x82, 0x40, 0x29, 0x13, 0x81, 0x41, 0xa2,
    ];
    let msg = [
        0x00, 0x1d, 0x0c, 0x23, 0x12, 0x87, 0xc1, 0x18, 0x27, 0x84, 0x55, 0x4c, 0xa3, 0xa2, 0x19,
        0x08,
    ];
    let aad = [];

    // ANCHOR: single_shot_encrypt
    let (ciphertext, tag) =
        aead::encrypt(Algorithm::Chacha20Poly1305, &key, &msg, &iv, &aad).unwrap();
    // ANCHOR_END: single_shot_encrypt

    // ANCHOR: single_shot_decrypt
    let msg_ = aead::decrypt(
        Algorithm::Chacha20Poly1305,
        &key,
        &ciphertext,
        &tag,
        &iv,
        &aad,
    )
    .unwrap();
    // ANCHOR_END: single_shot_decrypt

    assert_eq!(&msg[..], &msg_[..]);
    // ANCHOR_END: single_shot
}