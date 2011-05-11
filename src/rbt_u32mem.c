/*
 * Copyright 2010 Red Hat Inc., Durham, North Carolina.
 * All Rights Reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Authors:
 *      "Daniel Kopecek" <dkopecek@redhat.com>
 */
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "helpers.h"
#include "rbt_common.h"
#include "rbt_u32mem.h"

static struct rbt_node *rbt_u32mem_node_alloc(void)
{
        struct rbt_node *n;

        if (posix_memalign((void **)(void *)(&n), sizeof(void *),
                           sizeof (struct rbt_node) + sizeof (struct rbt_u32mem_node)) != 0)
        {
                abort ();
        }

        n->_chld[0] = NULL;
        n->_chld[1] = NULL;

        return(n);
}

static void rbt_u32mem_node_free(struct rbt_node *n)
{
        if (n != NULL)
                free(rbt_node_ptr(n));
}

rbt_t *rbt_u32mem_new (void)
{
        return rbt_new(RBT_U32MEMKEY);
}

static void rbt_u32mem_free_callback(struct rbt_u32mem_node *n)
{
        free(n->key);
        /* node memory is freed by rbt_free */
}

void rbt_u32mem_free (rbt_t *rbt)
{
        rbt_free(rbt, (void(*)(void *))&rbt_u32mem_free_callback);
}

void rbt_u32mem_free_cb (rbt_t *rbt, void (*callback)(struct rbt_u32mem_node *))
{
        rbt_free(rbt, (void(*)(void *))callback);
}

int rbt_u32mem_add(rbt_t *rbt, uint32_t *key, size_t keylen, void *ptr)
{
        struct rbt_node fake;
        register struct rbt_node *h[4];
        register uint8_t dvec;
        register uint32_t *n_key, *u_key;
        register int cmp;

        u_key = key;
        rbt_wlock(rbt);

        /*
         * Fake node
         */
        fake._chld[RBT_NODE_SL] = NULL;
        fake._chld[RBT_NODE_SR] = rbt->root;
        rbt_node_setcolor(&fake, RBT_NODE_CB);

        /*
         * Prepare node history stack & direction history vector
         */
        h[3] = NULL;
        h[2] = NULL;
        h[1] = &fake;
        h[0] = rbt->root;
        dvec = RBT_NODE_SR;

        for (;;) {
                if (rbt_node_ptr(h[0]) == NULL) {
                        /*
                         * Allocate new node aligned to sizeof(void *) bytes.
                         * On most systems, malloc already returns aligned
                         * memory but we want to ensure that its aligned using
                         * posix_memalign(3).
                         */
                        h[0] = rbt_u32mem_node_alloc ();

                        rbt_u32mem_node_key(h[0])    = key;
                        rbt_u32mem_node_keylen(h[0]) = keylen;
                        rbt_u32mem_node_ptr(h[0])    = ptr;

                        /*
                         * Set the new node as the child of the last node we've
                         * visited. The last direction is the lowest bit of the
                         * direction history vector.
                         */
                        rbt_node_setptr(rbt_node_ptr(h[1])->_chld[dvec & 1], h[0]);
                        rbt_node_setcolor(h[0], RBT_NODE_CR);

                        /*
                         * Since we are inserting a new red node, we need to fix
                         * a red violation if the parent of the new node is red.
                         */
                        if (rbt_node_getcolor(h[1]) == RBT_NODE_CR) {
                                rbt_redfix(h, dvec, rbt_node_ptr(h[3])->_chld[(dvec >> 2) & 1]);
                        }

                        /*
                         * Update root node and color it black
                         */
                        rbt->root = fake._chld[RBT_NODE_SR];
                        rbt_node_setcolor(rbt->root, RBT_NODE_CB);

                        /*
                         * Update node counter
                         */
                        ++(rbt->size);

                        break;
                } else if (rbt_node_getcolor(rbt_node_ptr(h[0])->_chld[0]) == RBT_NODE_CR &&
                           rbt_node_getcolor(rbt_node_ptr(h[0])->_chld[1]) == RBT_NODE_CR)
                {
                        /*
                         * Color switch
                         */
                        rbt_node_setcolor(h[0], RBT_NODE_CR);
                        rbt_node_setcolor(rbt_node_ptr(h[0])->_chld[0], RBT_NODE_CB);
                        rbt_node_setcolor(rbt_node_ptr(h[0])->_chld[1], RBT_NODE_CB);

                        /*
                         * Fix red violation
                         */
                        if (rbt_node_getcolor(h[1]) == RBT_NODE_CR) {
                                rbt_redfix(h, dvec, rbt_node_ptr(h[3])->_chld[(dvec >> 2) & 1]);
                        }
                } else if (rbt_node_getcolor(h[0]) == RBT_NODE_CR &&
                           rbt_node_getcolor(h[1]) == RBT_NODE_CR)
                {
                        /*
                         * Fix red violation
                         */
                        rbt_redfix(h, dvec, rbt_node_ptr(h[3])->_chld[(dvec >> 2)&1]);
                }

                n_key = rbt_u32mem_node_key(h[0]);
                cmp   = memcmp(u_key, n_key, sizeof(uint32_t) * rbt_u32mem_node_keylen(h[0]));

                if (cmp < 0) {
                        rbt_hpush4(h, rbt_node_ptr(h[0])->_chld[RBT_NODE_SL]);
                        dvec = (dvec << 1) | RBT_NODE_SL;
                } else if (cmp > 0) {
                        rbt_hpush4(h, rbt_node_ptr(h[0])->_chld[RBT_NODE_SR]);
                        dvec = (dvec << 1) | RBT_NODE_SR;
                } else {
                        /*
                         * Update root node and color it black
                         */
                        rbt->root = fake._chld[RBT_NODE_SR];
                        rbt_node_setcolor(rbt->root, RBT_NODE_CB);

                        rbt_wunlock(rbt);
                        return (1);
                }
        }

        rbt_wunlock(rbt);
        return (0);
}

int rbt_u32mem_get(rbt_t *rbt, uint32_t *key, size_t keylen, void **ptr)
{
        int r;
        struct rbt_u32mem_node *n = NULL;

        r = rbt_u32mem_getnode(rbt, key, keylen, &n);

        if (r == 0)
                *ptr = n->ptr;

        return(r);
}

int rbt_u32mem_getnode(rbt_t *rbt, uint32_t *key, size_t keylen, struct rbt_u32mem_node **node)
{
        int r;
        register struct rbt_node *n;
        register uint32_t *u_key, *n_key;
        register int cmp;

        u_key = (uint32_t *)key;
        r = -1;
        rbt_rlock(rbt);
        n = rbt_node_ptr(rbt->root);

        while (n != NULL) {
                assert(keylen == rbt_u32mem_node_keylen(n));
                n_key = rbt_u32mem_node_key(n);
                cmp   = memcmp(u_key, n_key, sizeof(uint32_t) * keylen);

                if (cmp < 0)
                        n = rbt_node_ptr(n->_chld[RBT_NODE_SL]);
                else if (cmp > 0)
                        n = rbt_node_ptr(n->_chld[RBT_NODE_SR]);
                else {
                        r = 0;
                        *node = (struct rbt_u32mem_node *)(n->_node);
                        break;
                }
        }

        rbt_runlock(rbt);
        return (r);
}
