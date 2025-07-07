#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pgapack.h"

#define N_ITEMS 1000

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

void print_node (rb_node_t *node)
{
    struct content *content = node->content;
    printf ("visit %d\n", content->i);
}

void check_consistency (rb_node_t *node)
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
    rb_walk (tree.root, count_black, check_consistency, NULL);
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
    for (i=0; i<N_ITEMS; i++) {
        scrambled [i] = i;
        rb_node_t *node = create_node (i);
        rb_node_t *parent = NULL;
        rb_node_t *found;
        dir_t dir = RB_LEFT;
        found = rb_search (&tree, node->content, &parent);
        assert (found == NULL);
        if (parent != NULL) {
            dir = cmp (node->content, parent->content) > 0;
        }
        rb_insert (&tree, node, parent, dir);
    }
    for (i=0; i<N_ITEMS; i++) {
        struct content c, *content;
        c.i = i;
        rb_node_t *parent = NULL;
        rb_node_t *n = rb_search (&tree, &c, &parent);
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
    rb_walk (tree.root, NULL, print_node, NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL);
    assert (tree.root->parent == NULL);
    for (i=0; i<N_ITEMS; i++) {
        struct content c;
        c.i = i;
        rb_node_t *parent = NULL;
        rb_node_t *n = rb_search (&tree, &c, &parent);
        printf ("n: %d\n", i);
        rb_walk (tree.root, count_black, check_consistency, NULL);
        assert (n != NULL);
        rb_remove (&tree, n);
        rb_walk (tree.root, count_black, check_consistency, NULL);
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
        rb_node_t *parent = NULL;
        rb_node_t *found;
        dir_t dir = RB_LEFT;
        found = rb_search (&tree, node->content, &parent);
        assert (found == NULL);
        if (parent != NULL) {
            dir = cmp (node->content, parent->content) > 0;
        }
        rb_insert (&tree, node, parent, dir);
    }
    /* Check consistency after insert */
    assert (tree.root->parent == NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL);
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
        rb_walk (tree.root, count_black, check_consistency, NULL);
    }
    /* Find nodes with one child and remove */
    left_count = right_count = 0;
    rb_walk (tree.root, NULL, NULL, rm_single_child_node);
    assert (left_count && right_count);
    rb_walk (tree.root, count_black, check_consistency, NULL);
    /* Finally loop over remaining nodes and remove again in random order */
    PGAShuffle (ctx, scrambled, N_ITEMS);
    for (i=0; i<N_ITEMS; i++) {
        int idx = scrambled [i];
        struct content c;
        rb_node_t *parent = NULL;
        rb_node_t *n;
        c.i = idx;
        n = rb_search (&tree, &c, &parent);
        rb_walk (tree.root, count_black, check_consistency, NULL);
        if (n != NULL) {
            printf ("Remove: %d\n", idx);
            rb_remove (&tree, n);
            rb_walk (tree.root, count_black, check_consistency, NULL);
            assert (n->content != NULL);
            free (n->content);
            memset (n, 0, sizeof (*n));
            free (n);
        }
    }
}
