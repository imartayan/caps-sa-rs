unsafe extern "C" {
    pub(crate) fn _build_sa_u8(text: *const u8, len: u32, ext_mem: bool, sa: *mut u32);
    pub(crate) fn _build_sa_u32(text: *const u32, len: u32, ext_mem: bool, sa: *mut u32);
    pub(crate) fn _build_sa_u64(text: *const u64, len: u32, ext_mem: bool, sa: *mut u32);

    pub(crate) fn _build_large_sa_u8(text: *const u8, len: u64, ext_mem: bool, sa: *mut u64);
    pub(crate) fn _build_large_sa_u32(text: *const u32, len: u64, ext_mem: bool, sa: *mut u64);
    pub(crate) fn _build_large_sa_u64(text: *const u64, len: u64, ext_mem: bool, sa: *mut u64);

    pub(crate) fn _build_sa_lcp_u8(
        text: *const u8,
        len: u32,
        ext_mem: bool,
        sa: *mut u32,
        lcp: *mut u32,
    );
    pub(crate) fn _build_sa_lcp_u32(
        text: *const u32,
        len: u32,
        ext_mem: bool,
        sa: *mut u32,
        lcp: *mut u32,
    );
    pub(crate) fn _build_sa_lcp_u64(
        text: *const u64,
        len: u32,
        ext_mem: bool,
        sa: *mut u32,
        lcp: *mut u32,
    );

    pub(crate) fn _build_large_sa_lcp_u8(
        text: *const u8,
        len: u64,
        ext_mem: bool,
        sa: *mut u64,
        lcp: *mut u64,
    );
    pub(crate) fn _build_large_sa_lcp_u32(
        text: *const u32,
        len: u64,
        ext_mem: bool,
        sa: *mut u64,
        lcp: *mut u64,
    );
    pub(crate) fn _build_large_sa_lcp_u64(
        text: *const u64,
        len: u64,
        ext_mem: bool,
        sa: *mut u64,
        lcp: *mut u64,
    );
}

macro_rules! build_sa_t {
    ($($t:ty),+) => {$(
        paste::paste! {
            #[doc = "Builds a suffix array from `text` encoded as `" $t "` and stores it in `sa`."]
            #[doc = ""]
            #[doc = "It uses external memory if `ext_mem` is set to `true`."]
            #[doc = "For texts larger than `2^32` characters, use [`build_large_sa_" $t "`] instead."]
            #[inline(always)]
            pub fn [<build_sa_$t>](text: &[$t], sa: &mut Vec<u32>, ext_mem: bool) {
                sa.reserve(text.len());
                unsafe {
                    [<_build_sa_$t>](text.as_ptr(), text.len() as u32, ext_mem, sa.as_mut_ptr());
                    sa.set_len(text.len());
                }
            }

            #[doc = "Builds a suffix array from large `text` encoded as `" $t "` and stores it in `sa`."]
            #[doc = ""]
            #[doc = "It uses external memory if `ext_mem` is set to `true`."]
            #[doc = "For texts smaller than `2^32` characters, [`build_sa_" $t "`] is recommended."]
            #[inline(always)]
            pub fn [<build_large_sa_$t>](text: &[$t], sa: &mut Vec<u64>, ext_mem: bool) {
                sa.reserve(text.len());
                unsafe {
                    [<_build_large_sa_$t>](text.as_ptr(), text.len() as u64, ext_mem, sa.as_mut_ptr());
                    sa.set_len(text.len());
                }
            }

            #[doc = "Builds a suffix array and LCP array from `text` encoded as `" $t "`, stores the suffix array in `sa` and the LCP array in `lcp`."]
            #[doc = ""]
            #[doc = "It uses external memory if `ext_mem` is set to `true`."]
            #[doc = "For texts larger than `2^32` characters, use [`build_large_sa_lcp_" $t "`] instead."]
            #[inline(always)]
            pub fn [<build_sa_lcp_$t>](text: &[$t], sa: &mut Vec<u32>, lcp: &mut Vec<u32>, ext_mem: bool) {
                if ext_mem {
                    unimplemented!();
                }
                sa.reserve(text.len());
                lcp.reserve(text.len());
                unsafe {
                    [<_build_sa_lcp_$t>](text.as_ptr(), text.len() as u32, ext_mem, sa.as_mut_ptr(), lcp.as_mut_ptr());
                    sa.set_len(text.len());
                    lcp.set_len(text.len());
                }
            }

            #[doc = "Builds a suffix array and LCP array from large `text` encoded as `" $t "`, stores the suffix array in `sa` and the LCP array in `lcp`."]
            #[doc = ""]
            #[doc = "It uses external memory if `ext_mem` is set to `true`."]
            #[doc = "For texts smaller than `2^32` characters, [`build_sa_lcp_" $t "`] is recommended."]
            #[inline(always)]
            pub fn [<build_large_sa_lcp_$t>](text: &[$t], sa: &mut Vec<u64>, lcp: &mut Vec<u64>, ext_mem: bool) {
                if ext_mem {
                    unimplemented!();
                }
                sa.reserve(text.len());
                lcp.reserve(text.len());
                unsafe {
                    [<_build_large_sa_lcp_$t>](text.as_ptr(), text.len() as u64, ext_mem, sa.as_mut_ptr(), lcp.as_mut_ptr());
                    sa.set_len(text.len());
                    lcp.set_len(text.len());
                }
            }
        }
    )*}
}

build_sa_t!(u8, u32, u64);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_sa_in_mem() {
        let text = [6, 7, 4, 5, 2, 3, 6, 7, 4, 5, 2, 3, 6, 7, 4, 5, 2, 3];
        let mut sa = Vec::new();
        build_sa_u32(&text, &mut sa, false);
        assert_eq!(
            sa,
            vec![16, 10, 4, 17, 11, 5, 14, 8, 2, 15, 9, 3, 12, 6, 0, 13, 7, 1]
        );

        let text: Vec<_> = (0..10_000).rev().collect();
        let mut sa = Vec::new();
        build_sa_u32(&text, &mut sa, false);
        assert_eq!(sa, text);
    }

    #[test]
    fn build_sa_ext_mem() {
        let text = [6, 7, 4, 5, 2, 3, 6, 7, 4, 5, 2, 3, 6, 7, 4, 5, 2, 3];
        let mut sa = Vec::new();
        build_sa_u32(&text, &mut sa, true);
        assert_eq!(
            sa,
            vec![16, 10, 4, 17, 11, 5, 14, 8, 2, 15, 9, 3, 12, 6, 0, 13, 7, 1]
        );

        let text: Vec<_> = (0..10_000).rev().collect();
        let mut sa = Vec::new();
        build_sa_u32(&text, &mut sa, true);
        for (i, x) in text.iter().copied().enumerate() {
            assert_eq!(sa[i], x);
        }
        assert_eq!(sa, text);
    }

    #[test]
    fn build_sa_lcp_in_mem() {
        let text = [6, 7, 4, 5, 2, 3, 6, 7, 4, 5, 2, 3, 6, 7, 4, 5, 2, 3];
        let mut sa = Vec::new();
        let mut lcp: Vec<u32> = Vec::new();
        build_sa_lcp_u32(&text, &mut sa, &mut lcp, false);
        assert_eq!(
            sa,
            vec![16, 10, 4, 17, 11, 5, 14, 8, 2, 15, 9, 3, 12, 6, 0, 13, 7, 1]
        );
        assert_eq!(
            lcp,
            vec![0, 2, 8, 0, 1, 7, 0, 4, 10, 0, 3, 9, 0, 6, 12, 0, 5, 11]
        );
    }
}
