#include "Suffix_Array.hpp"
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <string>

#ifdef __cplusplus
extern "C" {
#endif

void _build_sa_u8(const uint8_t *text, const uint32_t len, const bool ext_mem,
                  uint32_t *sa);
void _build_sa_u32(const uint32_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *sa);
void _build_sa_u64(const uint64_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *sa);

void _build_large_sa_u8(const uint8_t *text, const uint64_t len,
                        const bool ext_mem, uint64_t *sa);
void _build_large_sa_u32(const uint32_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *sa);
void _build_large_sa_u64(const uint64_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *sa);

void _build_sa_lcp_u8(const uint8_t *text, const uint32_t len,
                      const bool ext_mem, uint32_t *sa, uint32_t *lcp);
void _build_sa_lcp_u32(const uint32_t *text, const uint32_t len,
                       const bool ext_mem, uint32_t *sa, uint32_t *lcp);
void _build_sa_lcp_u64(const uint64_t *text, const uint32_t len,
                       const bool ext_mem, uint32_t *sa, uint32_t *lcp);

void _build_large_sa_lcp_u8(const uint8_t *text, const uint64_t len,
                            const bool ext_mem, uint64_t *sa, uint64_t *lcp);
void _build_large_sa_lcp_u32(const uint32_t *text, const uint64_t len,
                             const bool ext_mem, uint64_t *sa, uint64_t *lcp);
void _build_large_sa_lcp_u64(const uint64_t *text, const uint64_t len,
                             const bool ext_mem, uint64_t *sa, uint64_t *lcp);

#ifdef __cplusplus
}
#endif

template <typename T_seq_, typename T_idx_>
inline void _build_sa_lcp(const T_seq_ *text, const T_idx_ len,
                          const bool ext_mem, T_idx_ *sa,
                          T_idx_ *lcp = nullptr) {
  std::srand(std::time(0));
  const std::string ext_mem_path =
      "caps_sa_bucket_" + std::to_string(std::rand() % 512) + "_";
  const T_idx_ subproblem_count = 0;
  const T_idx_ max_context = 0;
  const bool output_lcp = lcp != nullptr;

  CaPS_SA::Suffix_Array<T_seq_, T_idx_> suf_arr(text, len, ext_mem,
                                                ext_mem_path, subproblem_count,
                                                max_context, output_lcp);
  ext_mem ? suf_arr.construct_ext_mem() : suf_arr.construct();

  std::memcpy(reinterpret_cast<void *>(sa), suf_arr.SA(), len * sizeof(T_idx_));
  if (output_lcp) {
    std::memcpy(reinterpret_cast<void *>(lcp), suf_arr.LCP(),
                len * sizeof(T_idx_));
  }
  if (ext_mem) {
    suf_arr.remove_extmem_partitions();
  }
}

void _build_sa_u8(const uint8_t *text, const uint32_t len, const bool ext_mem,
                  uint32_t *sa) {
  _build_sa_lcp(reinterpret_cast<const char *>(text), len, ext_mem, sa);
}

void _build_sa_u32(const uint32_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *sa) {
  _build_sa_lcp(text, len, ext_mem, sa);
}

void _build_sa_u64(const uint64_t *text, const uint32_t len, const bool ext_mem,
                   uint32_t *sa) {
  _build_sa_lcp(text, len, ext_mem, sa);
}

void _build_large_sa_u8(const uint8_t *text, const uint64_t len,
                        const bool ext_mem, uint64_t *sa) {
  _build_sa_lcp(reinterpret_cast<const char *>(text), len, ext_mem, sa);
}

void _build_large_sa_u32(const uint32_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *sa) {
  _build_sa_lcp(text, len, ext_mem, sa);
}

void _build_large_sa_u64(const uint64_t *text, const uint64_t len,
                         const bool ext_mem, uint64_t *sa) {
  _build_sa_lcp(text, len, ext_mem, sa);
}

void _build_sa_lcp_u8(const uint8_t *text, const uint32_t len,
                      const bool ext_mem, uint32_t *sa, uint32_t *lcp) {
  _build_sa_lcp(reinterpret_cast<const char *>(text), len, ext_mem, sa, lcp);
}

void _build_sa_lcp_u32(const uint32_t *text, const uint32_t len,
                       const bool ext_mem, uint32_t *sa, uint32_t *lcp) {
  _build_sa_lcp(text, len, ext_mem, sa, lcp);
}

void _build_sa_lcp_u64(const uint64_t *text, const uint32_t len,
                       const bool ext_mem, uint32_t *sa, uint32_t *lcp) {
  _build_sa_lcp(text, len, ext_mem, sa, lcp);
}

void _build_large_sa_lcp_u8(const uint8_t *text, const uint64_t len,
                            const bool ext_mem, uint64_t *sa, uint64_t *lcp) {
  _build_sa_lcp(reinterpret_cast<const char *>(text), len, ext_mem, sa, lcp);
}

void _build_large_sa_lcp_u32(const uint32_t *text, const uint64_t len,
                             const bool ext_mem, uint64_t *sa, uint64_t *lcp) {
  _build_sa_lcp(text, len, ext_mem, sa, lcp);
}

void _build_large_sa_lcp_u64(const uint64_t *text, const uint64_t len,
                             const bool ext_mem, uint64_t *sa, uint64_t *lcp) {
  _build_sa_lcp(text, len, ext_mem, sa, lcp);
}
