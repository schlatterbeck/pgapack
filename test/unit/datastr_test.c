#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pgapack.h"

#define N_ITEMS 1000
#define N_SMALL 10

static rb_tree_t tree;
static int scrambled [N_ITEMS];
static int left_count  = 0;
static int right_count = 0;

struct content {
    int i;
    int nblack;
};

int cmp (const void *a, const void *b)
{
    const struct content *v1 = a;
    const struct content *v2 = b;
    return v1->i < v2->i ? -1 : (v1->i > v2->i ? 1 : 0);
}

void print_node (rb_node_t *node, void *x)
{
    struct content *content = node->content;
    printf ("visit %d\n", content->i);
}

void check_consistency (rb_node_t *node, void *x)
{
    struct content *content = node->content;
    if (node->child [0] != NULL) {
        assert (node->child [0]->parent == node);
    }
    if (node->child [1] != NULL) {
        assert (node->child [1]->parent == node);
    }
    /* leaf */
    if (node->child [0] == NULL && node->child [1] == NULL) {
        rb_node_t *ll = rb_left_leaf (node);
        if (ll != NULL) {
            struct content *llcontent = ll->content;
            assert (llcontent->nblack == content->nblack);
        }
    }
}

void count_black (rb_node_t *node)
{
    struct content *content = node->content;
    int is_black = (node->color == RB_BLACK);
    if (node->parent == NULL) {
        content->nblack = is_black;
    } else {
        struct content *pcontent = node->parent->content;
        content->nblack = pcontent->nblack + is_black;
        if (node->parent->color == RB_RED) {
            assert (node->color == RB_BLACK);
        }
    }
}

void rm_single_child_node (rb_node_t *node)
{
    struct content *c;
    c = node->content;
    if (node->child [0] != NULL && node->child [1] == NULL) {
        left_count++;
        printf ("Remove: %d\n", c->i);
        rb_remove (&tree, node);
        free (c);
        /* We may not memset the node to zero during walk */
        free (node);
    }
    if (node->child [1] != NULL && node->child [0] == NULL) {
        right_count++;
        printf ("Remove: %d\n", c->i);
        rb_remove (&tree, node);
        free (c);
        /* We may not memset the node to zero during walk */
        free (node);
    }
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
}

rb_node_t *create_node (int i)
{
    rb_node_t *node = malloc (sizeof *node);
    struct content *content = malloc (sizeof (*content));

    if (node == NULL) {
        fprintf (stderr, "malloc node failed\n");
        exit (23);
    }
    memset (node, 0, sizeof (*node));
    if (content == NULL) {
        fprintf (stderr, "malloc content failed\n");
        exit (23);
    }
    memset (content, 0, sizeof (*content));
    content->i = i;
    content->nblack = -1;
    node->content = content;
    return node;
}

int main (int argc, char **argv)
{
    int i, j;
    PGAContext *ctx = NULL;
    tree.root = NULL;
    tree.cmp  = cmp;
    rb_node_t *fl = NULL;
    for (i=0; i<N_ITEMS; i++) {
        scrambled [i] = i;
        rb_node_t *node = create_node (i);
        rb_node_t *found;
        found = rb_search (&tree, node->content, NULL);
        assert (found == NULL);
        rb_insert (&tree, node);
    }
    i = 0;
    for (fl = rb_first (&tree), i=0; fl; fl = rb_next (fl), i++) {
        struct content *c = fl->content;
        assert (i == c->i);
    }
    assert (i == N_ITEMS);
    for (fl = rb_last (&tree), i=N_ITEMS - 1; fl; fl = rb_prev (fl), i--) {
        struct content *c = fl->content;
        assert (i == c->i);
    }
    assert (i == -1);
    for (i=0; i<N_ITEMS; i++) {
        struct content c, *content;
        c.i = i;
        rb_node_t *parent = NULL;
        rb_node_t *n = rb_search (&tree, &c, NULL);
        assert (n != NULL);
        content = n->content;
        assert (content->i == i);
        parent = n;
        for (j=0; j<22; j++) {
            parent = parent->parent;
            if (parent == NULL) {
                break;
            }
        }
        assert (j < 20);
    }
    rb_walk (tree.root, NULL, print_node, NULL, NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    assert (tree.root->parent == NULL);
    for (i=0; i<N_ITEMS; i++) {
        struct content c;
        c.i = i;
        rb_node_t *n = rb_search (&tree, &c, NULL);
        printf ("n: %d\n", i);
        rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
        assert (n != NULL);
        rb_remove (&tree, n);
        rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
        assert (n->content != NULL);
        free (n->content);
        memset (n, 0, sizeof (*n));
        free (n);
    }
    /* Now do this in random order */
    ctx = PGACreate (&argc, argv, PGA_DATATYPE_BINARY, 10, PGA_MINIMIZE);
    PGARandom01 (ctx, 42);
    PGAShuffle (ctx, scrambled, N_ITEMS);
    for (i=0; i<N_ITEMS; i++) {
        int idx = scrambled [i];
        rb_node_t *node = create_node (idx);
        rb_node_t *found;
        found = rb_search (&tree, node->content, NULL);
        assert (found == NULL);
        rb_insert (&tree, node);
    }
    /* Check consistency after insert */
    assert (tree.root->parent == NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    /* Test removing root node */
    for (i=0; i<5; i++) {
        struct content *c;
        rb_node_t *n = tree.root;
        c = n->content;
        printf ("Remove: %d\n", c->i);
        rb_remove (&tree, n);
        free (c);
        memset (n, 0, sizeof (*n));
        free (n);
        assert (tree.root->parent == NULL);
        rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    }
    /* Find nodes with one child and remove */
    left_count = right_count = 0;
    rb_walk (tree.root, NULL, NULL, rm_single_child_node, NULL);
    assert (left_count && right_count);
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    /* Finally loop over remaining nodes and remove again in random order */
    PGAShuffle (ctx, scrambled, N_ITEMS);
    for (i=0; i<N_ITEMS; i++) {
        int idx = scrambled [i];
        struct content c;
        rb_node_t *n;
        c.i = idx;
        n = rb_search (&tree, &c, NULL);
        rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
        if (n != NULL) {
            printf ("Remove: %d\n", idx);
            rb_remove (&tree, n);
            rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
            assert (n->content != NULL);
            free (n->content);
            memset (n, 0, sizeof (*n));
            free (n);
        }
    }
    /* empty tree to begin */
    assert (tree.root == NULL);
    /* Now insert three times */
    for (i=0; i<N_SMALL; i++) {
        scrambled [i] = i;
    }
    PGAShuffle (ctx, scrambled, N_SMALL);
    for (i=0; i<N_SMALL; i++) {
        int idx = scrambled [i];
        rb_node_t *node = create_node (idx);
        rb_node_t *found;
        found = rb_search (&tree, node->content, NULL);
        assert (found == NULL);
        rb_insert (&tree, node);
    }
    assert (tree.root->parent == NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    PGAShuffle (ctx, scrambled, N_SMALL);
    for (i=0; i<N_SMALL; i++) {
        int idx = scrambled [i];
        rb_node_t *node = create_node (idx);
        rb_node_t *found = rb_search (&tree, node->content, NULL);
        assert (found != NULL);
        rb_insert (&tree, node);
    }
    assert (tree.root->parent == NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    PGAShuffle (ctx, scrambled, N_SMALL);
    for (i=0; i<N_SMALL; i++) {
        int idx = scrambled [i];
        rb_node_t *node = create_node (idx);
        rb_node_t *found;
        found = rb_search (&tree, node->content, NULL);
        assert (found != NULL);
        rb_insert (&tree, node);
    }
    assert (tree.root->parent == NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    rb_walk (tree.root, NULL, print_node, NULL, NULL);
    PGAShuffle (ctx, scrambled, N_SMALL);
    for (i=0; i<N_SMALL * 3; i++) {
        rb_node_t *n;
        int idx = scrambled [i % N_SMALL];
        struct content c;
        c.i = idx;
        n = rb_search (&tree, &c, NULL);
        printf ("Remove: %d\n", idx);
        assert (n != NULL);
        rb_remove (&tree, n);
        rb_walk (tree.root, count_black, check_consistency, NULL, NULL);
    }
}
