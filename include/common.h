#ifndef __COMMON_H_
#define __COMMON_H_

enum class TimeOrderEnum { RK1, RK2, RK3 };
enum class BC { PERIODIC, NEUMANN, DIRICHLET, CUSTOM };
enum class PLTTYPE { ASCII, BINARY, BOTH, NONE };
enum class POISSONTYPE { MKL, CG, GS, BICGSTAB, ICPCG };

#endif __COMMON_H_