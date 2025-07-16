/*
COPYRIGHT

The following is a notice of limited availability of the code, and disclaimer
which must be included in the prologue of the code and in all source listings
of the code.

(C) COPYRIGHT 2025 Dr. Ralf Schlatterbeck Open Source Consulting

Permission is hereby granted to use, reproduce, prepare derivative works, and
to redistribute to others. This software was authored by:

Ralf Schlatterbeck
DISCLAIMER

This computer code material was prepared, in part, as an account of work
sponsored by an agency of the United States Government. Neither the United
States, nor the University of Chicago, nor any of their employees, makes any
warranty express or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not infringe
privately owned rights.
*/

/*!***************************************************************************
* \file
* This file contains implementations of data structures used in other parts
* of the code. Much of the red-black tree implementation was taken from
* the examples in Wikipedia, with some bugs removed.
* \authors Authors:
*          Ralf Schlatterbeck
*****************************************************************************/


#include <stddef.h>
#include <assert.h>
#include "pgapack.h"

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

#define DIRECTION(N) ((N)->parent == NULL \
                     ? RB_LEFT \
                     : ((N) == (N)->parent->child [1] ? RB_RIGHT : RB_LEFT) \
                     )

/*!****************************************************************************
    \brief Rotate a subtree in red/black tree
    \ingroup internal
    \param  tree  the tree
    \param  sub   subtree
    \param  dir   direction of rotation
    \return The new root after rotation

******************************************************************************/

static rb_node_t *rotate_subtree (rb_tree_t *tree, rb_node_t *sub, dir_t dir)
{
    rb_node_t *sub_parent = sub->parent;
    rb_node_t *new_root   = sub->child [1 - dir];
    rb_node_t *new_child  = new_root->child [dir];

    sub->child [1 - dir] = new_child;
    if (new_child != NULL) {
        new_child->parent = sub;
    }
    new_root->child [dir] = sub;
    new_root->parent      = sub_parent;
    sub->parent           = new_root;
    if (sub_parent != NULL) {
        sub_parent->child [sub == sub_parent->child [1]] = new_root;
    } else {
        tree->root = new_root;
    }
    return new_root;
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*!****************************************************************************
    \brief Compute smallest item in tree
    \ingroup internal

    \param  tree  the tree
    \return First node in tree or NULL if tree is empty

******************************************************************************/

rb_node_t *rb_first (const rb_tree_t *tree)
{
    rb_node_t *n = tree->root;
    if (n == NULL) {
        return NULL;
    }
    while (n->child [RB_LEFT] != NULL) {
        n = n->child [RB_LEFT];
    }
    return n;
}

/*!****************************************************************************
    \brief Internal insertion routine
    \ingroup internal

    \param  tree   the tree
    \param  node   the node to insert
    \param  parent the parent of the node
    \param  dir    direction from parent

    \rst

    Description
    -----------

    Assumes that node is the corrent point of insertion.

    \endrst

******************************************************************************/
static void rb_insert_internal
    (rb_tree_t *tree, rb_node_t *node, rb_node_t *parent, dir_t dir)
{
    rb_node_t *c = NULL;
    node->color  = RB_RED;
    node->child [0] = node->child [1] = node->parent = NULL;
    if (parent == NULL) {
        if (tree->root == NULL) {
            tree->root = node;
            return;
        } else {
            c = tree->root;
        }
    } else {
        c = parent->child [dir];
    }
    if (c != NULL) {
        rb_node_t *n = c;
        if (n->child [RB_LEFT] == NULL) {
            parent = n;
            dir = RB_LEFT;
        } else {
            n = n->child [RB_LEFT];
            while (n->child [RB_RIGHT]) {
                n = n->child [RB_RIGHT];
            }
            parent = n;
            dir = RB_RIGHT;
        }
    }
    assert (parent->child [dir] == NULL);
    parent->child [dir] = node;
    node->parent = parent;
    do {
        rb_node_t *grandparent = parent->parent;
        rb_node_t *uncle;
        if (parent->color == RB_BLACK) {
            return;
        }
        if (grandparent == NULL) {
            parent->color = RB_BLACK;
            return;
        }
        dir = DIRECTION (parent);
        uncle = grandparent->child [1 - dir];
        if (uncle == NULL || uncle->color == RB_BLACK) {
            if (node == parent->child [1 - dir]) {
                rotate_subtree (tree, parent, dir);
                node = parent;
                parent = grandparent->child [dir];
            }
            rotate_subtree (tree, grandparent, 1 - dir);
            parent->color = RB_BLACK;
            grandparent->color = RB_RED;
            return;
        }
        parent->color = RB_BLACK;
        uncle->color = RB_BLACK;
        grandparent->color = RB_RED;
        node = grandparent;
    } while (NULL != (parent = node->parent));
}

/*!****************************************************************************
    \brief Insert node into tree
    \ingroup internal

    \param  tree  the tree
    \param  node  the node to insert

******************************************************************************/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
void rb_insert (rb_tree_t *tree, rb_node_t *node)
{
    rb_node_t *parent = NULL;
    /* If assertions are turned off, 'found' is unused, see #pragma above */
    rb_node_t *found = rb_search (tree, node->content, &parent);
    assert (found == NULL || tree->cmp (found->content, node->content) == 0);
    dir_t dir = parent != NULL
              ? tree->cmp (node->content, parent->content) > 0
              : RB_LEFT
              ;
    rb_insert_internal (tree, node, parent, dir);
}
#pragma GCC diagnostic pop

/*!****************************************************************************
    \brief Compute largest item in tree
    \ingroup internal

    \param  tree  the tree
    \return The last node in the tree or NULL if tree is empty

******************************************************************************/
rb_node_t *rb_last (const rb_tree_t *tree)
{
    rb_node_t *n = tree->root;
    if (n == NULL) {
        return NULL;
    }
    while (n->child [RB_RIGHT] != NULL) {
        n = n->child [RB_RIGHT];
    }
    return n;
}

/*!****************************************************************************
    \brief Compute left leaf in tree from given node
    \ingroup internal

    \param  node  the start node
    \return The left leaf from the given node

******************************************************************************/
rb_node_t *rb_left_leaf (rb_node_t *node)
{
    rb_node_t *n = node;
    while (n->parent != NULL && DIRECTION (n) == RB_LEFT) {
        n = n->parent;
    }
    if (n->parent == NULL) {
        return NULL;
    }
    n = n->parent->child [0];
    if (n == NULL) {
        return NULL;
    }
    while (n->child [1]) {
        n = n->child [1];
    }
    return n;
}

/*!****************************************************************************
    \brief Compute next node from given node
    \ingroup internal

    \param  node  the start node
    \return The next node after node

******************************************************************************/
rb_node_t *rb_next (rb_node_t *node)
{
    rb_node_t *n = node;
    if (n->child [1]) {
        n = n->child [1];
        while (n->child [0]) {
            n = n->child [0];
        }
        return n;
    } else if (n->parent != NULL) {
        dir_t dir = DIRECTION (n);
        while (n->parent != NULL && dir == RB_RIGHT) {
            n = n->parent;
            dir = DIRECTION (n);
        }
        if (n->parent != NULL) {
            assert (dir == RB_LEFT);
            return n->parent;
        }
    }
    return NULL;
}

/*!****************************************************************************
    \brief Compute previous node from given node
    \ingroup internal

    \param  node  the start node
    \return The previous node before node

******************************************************************************/
rb_node_t *rb_prev (rb_node_t *node)
{
    rb_node_t *n = node;
    if (n->child [0]) {
        n = n->child [0];
        while (n->child [1]) {
            n = n->child [1];
        }
        return n;
    } else if (n->parent != NULL) {
        dir_t dir = DIRECTION (n);
        while (n->parent != NULL && dir == RB_LEFT) {
            n = n->parent;
            dir = DIRECTION (n);
        }
        if (n->parent != NULL) {
            assert (dir == RB_RIGHT);
            return n->parent;
        }
    }
    return NULL;
}

/*!****************************************************************************
    \brief Remove node in tree
    \ingroup internal

    \param  tree  the tree
    \param  node  the node to remove

******************************************************************************/
void rb_remove (rb_tree_t *tree, rb_node_t *node)
{
    rb_node_t *parent = node->parent;
    rb_node_t *sibling;
    rb_node_t *close_nephew;
    rb_node_t *distant_nephew;
    dir_t dir = DIRECTION (node);

    if (node->child [0] != NULL && node->child [1] != NULL) {
        /* Find leftmost child in right subtree */
        rb_node_t *c = node->child [1];
        while (c->child [0] != NULL) {
            c = c->child [0];
        }
        /* Removing the child may need rebalancing */
        rb_remove (tree, c);
        /* Now put this child into the node's place
         * Note that parent of node may have changed!
         * This also means dir can have changed!
         */
        if (node->parent != NULL) {
            node->parent->child [DIRECTION (node)] = c;
        } else {
            tree->root = c;
        }
        c->parent    = node->parent;
        c->child [0] = node->child [0];
        c->child [1] = node->child [1];
        c->color     = node->color;
        if (c->child [0] != NULL) {
            c->child [0]->parent = c;
        }
        if (c->child [1] != NULL) {
            c->child [1]->parent = c;
        }
        return;
    } else if (node->child [0] != NULL) {
        assert (node->color == RB_BLACK);
        assert (node->child [0]->color == RB_RED);
        if (parent == NULL) {
            tree->root = node->child [0];
        } else {
            parent->child [dir] = node->child [0];
        }
        node->child [0]->parent = parent;
        node->child [0]->color  = RB_BLACK;
        return;
    } else if (node->child [1] != NULL) {
        assert (node->color == RB_BLACK);
        assert (node->child [1]->color == RB_RED);
        if (parent == NULL) {
            tree->root = node->child [1];
        } else {
            parent->child [dir] = node->child [1];
        }
        node->child [1]->parent = parent;
        node->child [1]->color  = RB_BLACK;
        return;
    } else if (node->parent == NULL) {
        /* Tree is empty */
        tree->root = NULL;
        return;
    } else if (node->color == RB_RED) {
        parent->child [dir] = NULL;
        return;
    }

    /* This is the complicated case where the node has no children
     * and we need to re-balance
     */

    parent->child [dir] = NULL;

    do {
        sibling = parent->child [1 - dir];
        distant_nephew = sibling->child [1 - dir];
        close_nephew = sibling->child [dir];
        if (sibling->color == RB_RED) {
            // Case #3
            assert (close_nephew->color == RB_BLACK);
            assert (distant_nephew->color == RB_BLACK);
            rotate_subtree (tree, parent, dir);
            parent->color = RB_RED;
            sibling->color = RB_BLACK;
            sibling = close_nephew;

            distant_nephew = sibling->child [1 - dir];
            if (distant_nephew != NULL && distant_nephew->color == RB_RED) {
                goto case_6;
            }
            close_nephew = sibling->child [dir];
            if (close_nephew != NULL && close_nephew->color == RB_RED) {
                goto case_5;
            }

            // Case #4
            sibling->color = RB_RED;
            parent->color = RB_BLACK;
            return;
        }

        // Case #1
        if (parent == NULL) {
            return;
        }

        if (distant_nephew != NULL && distant_nephew->color == RB_RED) {
            goto case_6;
        }

        if (close_nephew != NULL && close_nephew->color == RB_RED) {
            goto case_5;
        }

        if (parent->color == RB_RED) {
            // Case #4
            sibling->color = RB_RED;
            parent->color = RB_BLACK;
            return;
        }

        // Case #2: all black
        sibling->color = RB_RED;
        node = parent;
        dir = DIRECTION (node);
    } while (NULL != (parent = node->parent));
    // case #1:
    return;

    case_5:

	rotate_subtree (tree, sibling, 1 - dir);
	sibling->color = RB_RED;
	close_nephew->color = RB_BLACK;
	distant_nephew = sibling;
	sibling = close_nephew;

    case_6:

	rotate_subtree (tree, parent, dir);
	sibling->color = parent->color;
	parent->color = RB_BLACK;
	distant_nephew->color = RB_BLACK;
	return;
}

/*!****************************************************************************
    \brief Search node in tree
    \ingroup internal

    \param  tree   the tree
    \param  item   Content to search for
    \param  parent Optional parent of found (or not found) node
    \return pointer to node found or NULL if no node is found

    \rst

    Description
    -----------

    When inserting the parent is returned for the insertion position
    even if the node is not found. Specify NULL for the parent if this
    feature is not needed.

    \endrst

******************************************************************************/
rb_node_t *rb_search
    (const rb_tree_t *tree, const void *item, rb_node_t **parent)
{
    rb_node_t *n = tree->root;
    if (parent != NULL) {
        *parent = NULL;
    }
    while (n != NULL) {
        int c = tree->cmp (item, n->content);
        switch (c) {
        case -1:
            if (parent != NULL) {
                *parent = n;
            }
            n = n->child [0];
            break;
        case  0:
            if (parent != NULL) {
                *parent = n->parent;
            }
            return n;
            break;
        case  1:
            if (parent != NULL) {
                *parent = n;
            }
            n = n->child [1];
            break;
        }
    }
    return n;
}

/*!****************************************************************************
    \brief Walk the tree
    \ingroup internal

    \param  tree      the node to start from
    \param  pre_func  Function to call before descending or NULL
    \param  func      Function before descending to right subtree or NULL
    \param  post_func Function to call before descending or NULL
    \param  payload   Payload to pass to func

    Description
    -----------

    This is a recursive method that descends into the left subtree, then
    into the right subtree. The function func is called in the middle
    and will get the nodes in tree sort order.

    \endrst

******************************************************************************/
void rb_walk
    ( rb_node_t *node
    , void (*pre_func)(rb_node_t *)
    , void (*func)(rb_node_t *, void *payload)
    , void (*post_func)(rb_node_t *)
    , void *payload
    )
{
    if (node == NULL) {
        return;
    }
    if (pre_func != NULL) {
        pre_func (node);
    }
    rb_walk (node->child [0], pre_func, func, post_func, payload);
    if (func != NULL) {
        func (node, payload);
    }
    rb_walk (node->child [1], pre_func, func, post_func, payload);
    if (post_func != NULL) {
        post_func (node);
    }
}
