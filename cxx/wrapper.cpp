#include "Suffix_Array.hpp"
#include <cstring>
#include <fstream>
#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void _build_sa_u8(const uint8_t *text, const uint32_t len, const bool ext_mem,
                  uint32_t *dest);
void _build_sa_u32(const uint32_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *dest);
void _build_sa_u64(const uint64_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *dest);

void _build_large_sa_u8(const uint8_t *text, const uint64_t len,
                        const bool ext_mem, uint64_t *dest);
void _build_large_sa_u32(const uint32_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *dest);
void _build_large_sa_u64(const uint64_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *dest);

#ifdef __cplusplus
}
#endif

template <typename T_seq_, typename T_idx_>
inline void _build_sa(const T_seq_ *text, const T_idx_ len, const bool ext_mem,
                      T_idx_ *dest) {
  const std::string ext_mem_path = "caps_sa_bucket";
  const T_idx_ subproblem_count = 0;
  const T_idx_ max_context = 0;
  const bool output_lcp = false;

  CaPS_SA::Suffix_Array<T_seq_, T_idx_> suf_arr(text, len, ext_mem,
                                                ext_mem_path, subproblem_count,
                                                max_context, output_lcp);
  ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();

  if (!ext_mem) {
    std::memcpy(reinterpret_cast<void *>(dest), suf_arr.SA(),
                len * sizeof(T_idx_));
  } else {
    for (T_idx_ p_id = 0; p_id < suf_arr.p(); ++p_id) {
      std::ifstream bucket(suf_arr.SA_bucket_file_path(p_id));
      bucket.seekg(0, std::ios::end);
      std::streamsize size = bucket.tellg();
      bucket.seekg(0, std::ios::beg);
      bucket.read(reinterpret_cast<char *>(dest), size);
      bucket.close();
    }
    suf_arr.remove_extmem_partitions();
  }
}

void _build_sa_u8(const uint8_t *text, const uint32_t len, const bool ext_mem,
                  uint32_t *dest) {
  _build_sa(reinterpret_cast<const char *>(text), len, ext_mem, dest);
}

void _build_sa_u32(const uint32_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *dest) {
  _build_sa(text, len, ext_mem, dest);
}

void _build_sa_u64(const uint64_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *dest) {
  _build_sa(text, len, ext_mem, dest);
}

void _build_large_sa_u8(const uint8_t *text, const uint64_t len,
                        const bool ext_mem, uint64_t *dest) {
  _build_sa(reinterpret_cast<const char *>(text), len, ext_mem, dest);
}

void _build_large_sa_u32(const uint32_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *dest) {
  _build_sa(text, len, ext_mem, dest);
}

void _build_large_sa_u64(const uint64_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *dest) {
  _build_sa(text, len, ext_mem, dest);
}
