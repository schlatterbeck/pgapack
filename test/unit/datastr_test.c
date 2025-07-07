#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pgapack.h"

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

int main ()
{
    int i, j;
    rb_tree_t tree;
    tree.root = NULL;
    tree.cmp  = cmp;
    for (i=0; i<1000; i++) {
        rb_node_t *node = malloc (sizeof *node);
        struct content *content = malloc (sizeof (*content));
        rb_node_t *parent = NULL;
        rb_node_t *found;
        dir_t dir = RB_LEFT;
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
        found = rb_search (&tree, content, &parent);
        assert (found == NULL);
        if (parent != NULL) {
            dir = cmp (node->content, parent->content) > 0;
        }
        rb_insert (&tree, node, parent, dir);
    }
    for (i=0; i<1000; i++) {
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
        //printf ("depth: %d\n", j);
        assert (j < 20);
    }
    rb_walk (tree.root, NULL, print_node, NULL);
    rb_walk (tree.root, count_black, check_consistency, NULL);
    for (i=0; i<1000; i++) {
        struct content c;
        c.i = i;
        rb_node_t *parent = NULL;
        rb_node_t *n = rb_search (&tree, &c, &parent);
        printf ("n: %d\n", i);
        rb_walk (tree.root, count_black, check_consistency, NULL);
        assert (n != NULL);
        rb_remove (&tree, n);
        rb_walk (tree.root, count_black, check_consistency, NULL);
        memset (n, 0, sizeof (*n));
        free (n);
    }
}
