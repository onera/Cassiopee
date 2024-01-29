#include "Proto.h"

Tree::Tree(E_Int nelem) :
  nelem_(nelem),
  parent_(NULL), children_(NULL),
  level_(NULL), type_(NULL), state_(NULL)
{
  parent_ = (E_Int *)XMALLOC(nelem * sizeof(E_Int));
  children_ = (Children **)XMALLOC(nelem * sizeof(Children *));
  level_ = (E_Int *)XMALLOC(nelem * sizeof(E_Int));
  type_ = (E_Int *)XMALLOC(nelem * sizeof(E_Int));
  state_ = (E_Int *)XMALLOC(nelem * sizeof(E_Int));

  for (E_Int i = 0; i < nelem; i++) {
    parent_[i] = i;
    children_[i] = NULL;
    level_[i] = 0;
    type_[i] = TYPE_UNKNOWN;
    state_[i] = UNTOUCHED;
  }
}

void Tree::resize(E_Int new_elem)
{
  parent_ = (E_Int *)XRESIZE(parent_, new_elem * sizeof(E_Int));
  children_ = (Children **)XRESIZE(children_, new_elem * sizeof(Children *));
  level_ = (E_Int *)XRESIZE(level_, new_elem * sizeof(E_Int));
  type_ = (E_Int *)XRESIZE(type_, new_elem * sizeof(E_Int));
  state_ = (E_Int *)XRESIZE(state_, new_elem * sizeof(E_Int));

  for (E_Int i = nelem_; i < new_elem; i++) {
    parent_[i] = i;
    children_[i] = NULL;
    level_[i] = 0;
    type_[i] = TYPE_UNKNOWN;
    state_[i] = UNTOUCHED;
  }

  nelem_ = new_elem;
}

static
void print_children_node(Children *children)
{
  printf("[ ");
  for (E_Int i = 0; i < children->n; i++)
    printf("%d ", children->pc[i]);
  printf("]");
}

void Tree::print_elem_histo(E_Int elem)
{
  Children *children = children_[elem];

  while (children) {
    print_children_node(children);
    printf(" -> ");
    children = children->next;
  }
  puts("");
}

void Tree::print_children()
{
  for (E_Int i = 0; i < nelem_; i++) {
    if (!children_[i]) {
      //printf("---");
    } else {
      printf("%d -> ", i);
      Children *cur = children_[i];
      while (cur != NULL) {
        print_children_node(cur);
        if (cur->next)
          printf(" -> ");
        cur = cur->next;
      }
      puts("");
    }
  }
}

void Tree::print_face_types()
{
  for (E_Int i = 0; i < nelem_; i++) {
    printf("%d -> %d (%s) \n", i, type_[i], face_type_to_string(type_[i]));
  }
}

void Tree::print_cell_types()
{
  for (E_Int i = 0; i < nelem_; i++) {
    printf("%d -> %d (%s) \n", i, type_[i], cell_type_to_string(type_[i]));
  }
}

void Tree::set_new_elem(E_Int elem, E_Int type, E_Int level)
{
  assert(children_[elem] == NULL);
  parent_[elem] = elem;
  type_[elem] = type;
  level_[elem] = level;
  state_[elem] = UNTOUCHED;
}

void Tree::set_parent_elem(E_Int elem, E_Int nchildren, E_Int pos)
{
  Children *children = (Children *)XMALLOC(sizeof(Children));
  children->n = nchildren;
  children->pc[0] = elem;

  for (E_Int i = 1; i < nchildren; i++)
    children->pc[i] = pos + i - 1;
 
  children->next = children_[elem];
  children_[elem] = children;

  level_[elem] += 1;
  state_[elem] = REFINED;
}

void Tree::set_child_elem(E_Int cpos, E_Int parent, E_Int type, E_Int level,
  E_Int pos)
{
  E_Int id = cpos + pos - 1;
  assert(children_[id] == NULL);
  parent_[id] = parent;
  type_[id] = type;
  level_[id] = level;
  state_[id] = UNTOUCHED;
}

void Tree::compress(const std::vector<E_Int> &new_ids, E_Int new_nelem)
{
  E_Int *new_parent = (E_Int *)XMALLOC(new_nelem * sizeof(E_Int));
  Children **new_children = (Children **)XMALLOC(new_nelem * sizeof(Children *));
  E_Int *new_level = (E_Int *)XMALLOC(new_nelem * sizeof(E_Int));
  E_Int *new_type = (E_Int *)XMALLOC(new_nelem * sizeof(E_Int));
  E_Int *new_state = (E_Int *)XMALLOC(new_nelem * sizeof(E_Int));

  assert(nelem_ == (E_Int)new_ids.size());

  // Shift
  E_Int j = 0;
  for (E_Int i = 0; i < nelem_; i++) {
    if (state(i) != GONE) {
      new_parent[j] = parent_[i];
      new_children[j] = children_[i];
      new_level[j] = level_[i];
      new_type[j] = type_[i];
      new_state[j] = state_[i];
      j++;
    }
  }
  assert(j == new_nelem);

  // Renumber
  for (E_Int i = 0; i < new_nelem; i++) {
    new_parent[i] = new_ids[new_parent[i]];
    assert(new_parent[i] != -1);
    Children *cur = new_children[i];
    while (cur != NULL) {
      for (E_Int j = 0; j < cur->n; j++) {
        cur->pc[j] = new_ids[cur->pc[j]];
        assert(cur->pc[j] != -1);
      }

      cur = cur->next;
    }
  }

  XFREE(parent_);
  XFREE(level_);
  XFREE(type_);
  XFREE(state_);
  XFREE(children_);

  nelem_ = new_nelem;
  parent_ = new_parent;
  level_ = new_level;
  type_ = new_type;
  state_ = new_state;
  children_ = new_children;
}

void Tree::drop()
{
  XFREE(parent_);
  for (E_Int i = 0; i < nelem_; i++) {
    Children *cur = children_[i];

    while (cur) {
      Children *tmp = cur->next;
      XFREE(cur);
      cur = tmp;
    }
  }
  XFREE(children_);
  XFREE(state_);
  XFREE(level_);
  XFREE(type_);
}

void Tree::renumber(const std::vector<E_Int> &new_ids)
{
  E_Int *parent = (E_Int *)XMALLOC(nelem_ * sizeof(E_Int));
  E_Int *level = (E_Int *)XMALLOC(nelem_ * sizeof(E_Int));
  E_Int *type = (E_Int *)XMALLOC(nelem_ * sizeof(E_Int));
  E_Int *state = (E_Int *)XMALLOC(nelem_ * sizeof(E_Int));
  Children **children = (Children **)XMALLOC(nelem_ * sizeof(Children *));

  for (E_Int i = 0; i < nelem_; i++) {
    parent[new_ids[i]] = parent_[i];
    level[new_ids[i]] = level_[i];
    type[new_ids[i]] = type_[i];
    state[new_ids[i]] = state_[i];
    children[new_ids[i]] = children_[i];
  }

  for (E_Int i = 0; i < nelem_; i++) {
    parent[i] = new_ids[parent[i]];

    Children *cur = children[i];
    while (cur) {
      for (E_Int j = 0; j < cur->n; j++) {
        cur->pc[j] = new_ids[cur->pc[j]];
      }

      cur = cur->next;
    }
  }

  XFREE(parent_);
  XFREE(level_);
  XFREE(state_);
  XFREE(type_);
  XFREE(children_);

  parent_ = parent;
  level_ = level;
  state_ = state;
  type_ = type;
  children_ = children;
}

int Tree::check_numbering(E_Int N)
{
  E_Int *velems = (E_Int *)XCALLOC(nelem_, sizeof(E_Int));
  //for (E_Int i = 0; i < N; i++) velems[i] = 1;

  for (E_Int i = 0; i < N; i++) {
    velems[i] = 1;
    Children *C = children_[i];

    while (C) {
      for (E_Int j = 0; j < C->n; j++)
        velems[C->pc[j]] = 1;
      C = C->next;
    }
  }

  E_Int ok = 1;
  for (E_Int i = 0; i < N; i++) {
    if (velems[i] != 1) {
      ok = 0;
      break;
    }
  }

  XFREE(velems);

  return ok;
}
