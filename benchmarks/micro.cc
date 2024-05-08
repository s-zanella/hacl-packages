#include <inttypes.h>
#include <string.h>
#include "util.h"


static const int32_t libcrux_kyber_ntt_ZETAS_TIMES_MONTGOMERY_R[128U] = {
  (int32_t)-1044, (int32_t)-758,  (int32_t)-359,  (int32_t)-1517,
  (int32_t)1493,  (int32_t)1422,  (int32_t)287,   (int32_t)202,
  (int32_t)-171,  (int32_t)622,   (int32_t)1577,  (int32_t)182,
  (int32_t)962,   (int32_t)-1202, (int32_t)-1474, (int32_t)1468,
  (int32_t)573,   (int32_t)-1325, (int32_t)264,   (int32_t)383,
  (int32_t)-829,  (int32_t)1458,  (int32_t)-1602, (int32_t)-130,
  (int32_t)-681,  (int32_t)1017,  (int32_t)732,   (int32_t)608,
  (int32_t)-1542, (int32_t)411,   (int32_t)-205,  (int32_t)-1571,
  (int32_t)1223,  (int32_t)652,   (int32_t)-552,  (int32_t)1015,
  (int32_t)-1293, (int32_t)1491,  (int32_t)-282,  (int32_t)-1544,
  (int32_t)516,   (int32_t)-8,    (int32_t)-320,  (int32_t)-666,
  (int32_t)-1618, (int32_t)-1162, (int32_t)126,   (int32_t)1469,
  (int32_t)-853,  (int32_t)-90,   (int32_t)-271,  (int32_t)830,
  (int32_t)107,   (int32_t)-1421, (int32_t)-247,  (int32_t)-951,
  (int32_t)-398,  (int32_t)961,   (int32_t)-1508, (int32_t)-725,
  (int32_t)448,   (int32_t)-1065, (int32_t)677,   (int32_t)-1275,
  (int32_t)-1103, (int32_t)430,   (int32_t)555,   (int32_t)843,
  (int32_t)-1251, (int32_t)871,   (int32_t)1550,  (int32_t)105,
  (int32_t)422,   (int32_t)587,   (int32_t)177,   (int32_t)-235,
  (int32_t)-291,  (int32_t)-460,  (int32_t)1574,  (int32_t)1653,
  (int32_t)-246,  (int32_t)778,   (int32_t)1159,  (int32_t)-147,
  (int32_t)-777,  (int32_t)1483,  (int32_t)-602,  (int32_t)1119,
  (int32_t)-1590, (int32_t)644,   (int32_t)-872,  (int32_t)349,
  (int32_t)418,   (int32_t)329,   (int32_t)-156,  (int32_t)-75,
  (int32_t)817,   (int32_t)1097,  (int32_t)603,   (int32_t)610,
  (int32_t)1322,  (int32_t)-1285, (int32_t)-1465, (int32_t)384,
  (int32_t)-1215, (int32_t)-136,  (int32_t)1218,  (int32_t)-1335,
  (int32_t)-874,  (int32_t)220,   (int32_t)-1187, (int32_t)-1659,
  (int32_t)-1185, (int32_t)-1530, (int32_t)-1278, (int32_t)794,
  (int32_t)-1510, (int32_t)-854,  (int32_t)-870,  (int32_t)478,
  (int32_t)-108,  (int32_t)-308,  (int32_t)996,   (int32_t)991,
  (int32_t)958,   (int32_t)-1460, (int32_t)1522,  (int32_t)1628
};

static inline uint32_t
libcrux_kyber_arithmetic_get_n_least_significant_bits(uint8_t n, uint32_t value)
{
  return value & ((1U << (uint32_t)n) - 1U);
}

static inline int32_t
libcrux_kyber_arithmetic_montgomery_reduce(int32_t value)
{
  uint32_t t = libcrux_kyber_arithmetic_get_n_least_significant_bits(
                 (16U), (uint32_t)value) *
               (62209U);
  int16_t k =
    (int16_t)libcrux_kyber_arithmetic_get_n_least_significant_bits((16U), t);
  int32_t k_times_modulus = (int32_t)k * ((int32_t)3329);
  int32_t c = k_times_modulus >> (uint32_t)(16U);
  int32_t value_high = value >> (uint32_t)(16U);
  return value_high - c;
}

static inline int32_t
libcrux_kyber_arithmetic_montgomery_multiply_fe_by_fer(int32_t fe, int32_t fer)
{
  return libcrux_kyber_arithmetic_montgomery_reduce(fe * fer);
}

static inline void
libcrux_kyber_ntt_ntt_at_layer(size_t* zeta_i,
                               int32_t re[256U],
                               size_t layer,
                               size_t _initial_coefficient_bound,
                               int32_t ret[256U])
{
  // unsigned __int64 _start = __rdtsc();
  size_t step = (size_t)1U << (uint32_t)layer;
  for (size_t i0 = (size_t)0U; i0 < (size_t)128U >> (uint32_t)layer; i0++) {
    size_t round = i0;
    zeta_i[0U] = zeta_i[0U] + (size_t)1U;
    size_t offset = round * step * (size_t)2U;
    for (size_t i = offset; i < offset + step; i++) {
      size_t j = i;
      int32_t t = libcrux_kyber_arithmetic_montgomery_multiply_fe_by_fer(
        re[j + step], libcrux_kyber_ntt_ZETAS_TIMES_MONTGOMERY_R[zeta_i[0U]]);
      re[j + step] = re[j] - t;
      re[j] = re[j] + t;
    }
  }
  memcpy(ret, re, (size_t)256U * sizeof(int32_t));
  //_total += __rdtsc() - _start;
  //_count++;
}

static inline void
libcrux_kyber_ntt_ntt_at_layer_3(size_t* zeta_i,
                                 int32_t re[256U],
                                 size_t layer,
                                 int32_t ret[256U])
{
  int32_t ret0[256U];
  libcrux_kyber_ntt_ntt_at_layer(zeta_i, re, layer, (size_t)3U, ret0);
  memcpy(ret, ret0, (size_t)256U * sizeof(int32_t));
}

static inline int64_t
core_convert_num___core__convert__From_i32__for_i64__59__from(int32_t x)
{
  return x;
}

static inline int32_t
libcrux_kyber_arithmetic_barrett_reduce(int32_t value)
{
  int64_t t =
    core_convert_num___core__convert__From_i32__for_i64__59__from(value) *
      ((int64_t)20159) +
    (((int64_t)1 << (uint32_t)((int64_t)26)) >> 1U);
  int32_t quotient = (int32_t)(t >> (uint32_t)((int64_t)26));
  return value - quotient * ((int32_t)3329);
}

static inline void
libcrux_kyber_ntt_ntt_binomially_sampled_ring_element(int32_t re[256U],
                                                      int32_t ret[256U])
{
  size_t zeta_i = (size_t)1U;
  for (size_t i = (size_t)0U; i < (size_t)128U; i++) {
    size_t j = i;
    int32_t t = re[j + (size_t)128U] * (int32_t)-1600;
    re[j + (size_t)128U] = re[j] - t;
    re[j] = re[j] + t;
  }
  libcrux_kyber_ntt_ntt_at_layer_3(&zeta_i, re, (size_t)6U, re);
  libcrux_kyber_ntt_ntt_at_layer_3(&zeta_i, re, (size_t)5U, re);
  libcrux_kyber_ntt_ntt_at_layer_3(&zeta_i, re, (size_t)4U, re);
  libcrux_kyber_ntt_ntt_at_layer_3(&zeta_i, re, (size_t)3U, re);
  libcrux_kyber_ntt_ntt_at_layer_3(&zeta_i, re, (size_t)2U, re);
  libcrux_kyber_ntt_ntt_at_layer_3(&zeta_i, re, (size_t)1U, re);
  for (size_t i = (size_t)0U; i < ((size_t)256U); i++) {
    size_t i0 = i;
    int32_t uu____0 = libcrux_kyber_arithmetic_barrett_reduce(re[i0]);
    re[i0] = uu____0;
  }
  memcpy(ret, re, (size_t)256U * sizeof(int32_t));
}


static void
micro(benchmark::State& state)
{
  int32_t re[256U];
  int32_t ret[256U];

  generate_random((uint8_t *)re, (sizeof re)/(sizeof re[0]));
  
  for (auto _ : state) {
    libcrux_kyber_ntt_ntt_binomially_sampled_ring_element(re, ret);
  }
}

BENCHMARK(micro)->Setup(DoSetup);

BENCHMARK_MAIN();