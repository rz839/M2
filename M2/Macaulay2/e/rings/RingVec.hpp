#include "Ring.h"
#include "text-io.hpp"
#include <vector>
//#include "matrix.hpp"
//#include "geovec.hpp"
//#include "ringmap.hpp"

class Matrix;
class RingMap;

//#include "poly.hpp"
//  Notes: ring_elem's are treated as immutable objects: they are not changed,
//  and
// the fact that one cannot change is used throughout.


template <typename D>
void RingBase<D>::remove_vec_node(vec n) const
{
  deleteitem(n);
}

template <typename D>
void RingBase<D>::vec_sort(vecterm *&f) const
{
  // Internal routine to place the elements back in order after
  // an operation such as subvector.
  // Divide f into two lists of equal length, sort each,
  // then add them together.  This allows the same monomial
  // to appear more than once in 'f'.

  if (f == nullptr || f->next == nullptr) return;
  vecterm *f1 = nullptr;
  vecterm *f2 = nullptr;
  while (f != nullptr)
  {
    vecterm *t = f;
    f = f->next;
    t->next = f1;
    f1 = t;

    if (f == nullptr) break;
    t = f;
    f = f->next;
    t->next = f2;
    f2 = t;
  }

  vec_sort(f1);
  vec_sort(f2);
  add_vec_to(f1, f2);
  f = f1;
}

template <typename D>
int RingBase<D>::compare_vecs(vec v, vec w) const
{
  for (;; v = v->next, w = w->next)
  {
    if (v == nullptr)
    {
      if (w == nullptr) return 0;
      return -1;
    }
    if (w == nullptr) return 1;
    int cmp = v->comp - w->comp;
    if (cmp > 0) return cmp;
    if (cmp < 0) return cmp;
    cmp = this->compare_elems(v->coeff, w->coeff);
    if (cmp > 0) return cmp;
    if (cmp < 0) return cmp;
  }
}

template <typename D>
vec RingBase<D>::e_sub_i(int i) const
{
  ring_elem a = crtp()->from_long(1);
  return make_vec(i, a);
}

template <typename D>
vec RingBase<D>::make_vec(int r, ring_elem a) const
{
  if (crtp()->is_zero(a)) return nullptr;
  vec result = new_vec();
  result->next = 0;
  result->comp = r;
  result->coeff = a;
  return result;
}

template <typename D>
vec RingBase<D>::make_vec_from_array(int len, Nterm **array) const
{
  vec result = 0;
  for (int i = 0; i < len; i++)
  {
    if (array[i] != 0)
    {
      vec v = make_vec(i, array[i]);
      v->next = result;
      result = v;
    }
  }
  return result;
}

template <typename D>
vec RingBase<D>::copy_vec(const vecterm *v) const
{
  vecterm head;
  vec result = &head;
  for (const vecterm *p = v; p != 0; p = p->next)
  {
    vec w = new_vec();
    result->next = w;
    result = w;
    w->comp = p->comp;
    w->coeff = p->coeff;  // copy is not done
  }
  result->next = 0;
  return head.next;
}

template <typename D>
void RingBase<D>::remove_vec(vec v) const
{
  while (v != 0)
  {
    vec tmp = v;
    v = v->next;
    remove_vec_node(tmp);
  }
}

template <typename D>
void RingBase<D>::add_vec_to(vec &v, vec &w) const
{
  if (w == NULL) return;
  if (v == NULL)
    {
      v = w;
      w = NULL;
      return;
    }
  vecterm head;
  vec result = &head;
  while (true)
    if (v->comp < w->comp)
      {
        result->next = w;
        result = result->next;
        w = w->next;
        if (w == NULL)
          {
            result->next = v;
            v = head.next;
            return;
          }
      }
    else if (v->comp > w->comp)
      {
        result->next = v;
        result = result->next;
        v = v->next;
        if (v == NULL)
          {
            result->next = w;
            v = head.next;
            w = NULL;
            return;
          }
      }
    else
      {
        vec tmv = v;
        vec tmw = w;
        v = v->next;
        w = w->next;
        tmv->coeff = this->add(tmv->coeff, tmw->coeff);
        if (this->is_zero(tmv->coeff))
          {
            remove_vec_node(tmv);
          }
        else
          {
            result->next = tmv;
            result = result->next;
          }
        remove_vec_node(tmw);
        if (w == NULL)
          {
            result->next = v;
            v = head.next;
            return;
          }
        if (v == NULL)
          {
            result->next = w;
            v = head.next;
            w = NULL;
            return;
          }
      }
}

// Local Variables:
// compile-command: "make -C $M2BUILDDIR/Macaulay2/e "
// indent-tabs-mode: nil
// End:
