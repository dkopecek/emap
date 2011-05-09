bin_PROGRAMS= emap

emap_SOURCES=	emap.c		\
		emap.h		\
		pointdb.c	\
		pointdb.h	\
		helpers.c	\
		helpers.h	\
		rbt_common.c	\
		rbt_common.h	\
		rbt_u32mem.c	\
		rbt_u32mem.h

emap_CFLAGS= -DEMAP_POINT_DOUBLE -DRBT_IMPLITIC_LOCKING -fopenmp
emap_LDFLAGS= -lm
