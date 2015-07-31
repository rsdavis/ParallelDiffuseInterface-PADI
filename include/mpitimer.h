
#ifndef MPITIMER_H
#define MPITIMER_H

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

typedef struct {
    double start;
    double elapsed;
    int on;
} mpitimer;

static inline mpitimer * mpitimer_new()
{
    mpitimer * t;
    t = (mpitimer *) malloc(sizeof(mpitimer));
    t->on = 0;
    t->start = 0;
    t->elapsed = 0;
    return t;
}

static inline void mpitimer_delete(mpitimer * t)
{
    free(t);
}

static inline void mpitimer_start(mpitimer * const t)
{
    t->on = 1;
    t->start = MPI_Wtime();
}

static inline void mpitimer_stop(mpitimer * const t)
{
    if (t->on) t->elapsed += MPI_Wtime() - t->start;
    t->on = 0;
}

static inline int mpitimer_is_running(mpitimer * const t)
{
    return t->on;
}

static inline void mpitimer_reset(mpitimer * const t)
{
    t->on = 0;
    t->start = 0;
    t->elapsed = 0;
}

static inline double mpitimer_get_time(const mpitimer * const t)
{
    double time;
    if (t->on == 0) time = t->elapsed;
    else time = t->elapsed + MPI_Wtime() - t->start;
    return time; 
}

#ifdef __cplusplus
}
#endif

#endif
