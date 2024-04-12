/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "proto.h"

PyObject *K_XCORE::sweep(PyObject *self, PyObject *args)
{
  PyObject *MASTER, *SLAVE;

  if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVE)) {
    RAISE("Bad input.");
    return NULL;
  }

  // Master mesh
  E_Int ret, ni, nj, nk;
  K_FLD::FldArrayF *fm;
  K_FLD::FldArrayI *cnm;
  char *varString, *eltType;
  ret = K_ARRAY::getFromArray3(MASTER, varString, fm, ni, nj, nk, cnm, eltType);

  assert(ret == 2);

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  assert(posx != -1 && posy != -1 && posz != -1);

  posx++; posy++; posz++;

  E_Float *Xm = fm->begin(posx);
  E_Float *Ym = fm->begin(posy);
  E_Float *Zm = fm->begin(posz);
  E_Int npm = fm->getSize();

  // Slave mesh
  K_FLD::FldArrayF *fs;
  K_FLD::FldArrayI *cns;
  ret = K_ARRAY::getFromArray3(SLAVE, varString, fs, ni, nj, nk, cns, eltType);

  assert(ret == 2);

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);

  assert(posx != -1 && posy != -1 && posz != -1);

  posx++; posy++; posz++;

  E_Float *Xs = fs->begin(posx);
  E_Float *Ys = fs->begin(posy);
  E_Float *Zs = fs->begin(posz);
  E_Int nps = fs->getSize();

  // Init meshes
  Mesh *M = mesh_init(*cnm, Xm, Ym, Zm, npm);
  Mesh *S = mesh_init(*cns, Xs, Ys, Zs, nps);

  printf("Master faces:  " SF_D_ "\n", M->nc);
  printf("Master edges:  " SF_D_ "\n", M->nf);
  printf("Master points: " SF_D_ "\n", M->np);
  printf("Slave faces:   " SF_D_ "\n", S->nc);
  printf("Slave edges:   " SF_D_ "\n", S->nf);
  printf("Slave points:  " SF_D_ "\n", S->np);

  // Make points
  std::vector<point *> points;

  for (E_Int i = 0; i < M->np; i++) {
    point *p = new point(M->x[i], M->y[i], i, -1);
    points.push_back(p);
  }

  for (E_Int i = 0; i < S->np; i++) {
    point *p = new point(S->x[i], S->y[i], i+M->np, -1);
    points.push_back(p);
  }

  for (point *p : points)
    printf("%zu: %f %f\n", p->id, p->x, p->y);

  // Make segments
  std::vector<segment *> segs;

  for (E_Int i = 0; i < M->nf; i++) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(i, np, M);
    assert(np == 2);
    segment *s = new segment(points[pn[0]], points[pn[1]], i, i);
    if (!segment_is_lexico(s))
      std::swap(s->p, s->q);
    segs.push_back(s);
  }

  for (E_Int i = 0; i < S->nf; i++) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(i, np, S);
    assert(np == 2);
    segment *s = new segment(points[pn[0]+M->np], points[pn[1]+M->np], M->nf+i,
      M->nf+i);
    if (!segment_is_lexico(s))
      std::swap(s->p, s->q);
    segs.push_back(s);
  }

  for (segment *s: segs)
    printf("S%zu: P%zu P%zu\n", s->id, s->p->id, s->q->id);

  std::unordered_map<point *, std::vector<segment *>> Xsegs;

  _sweep(segs, points, Xsegs);

  return Py_None;
}


static
int _quicksort_partition(std::vector<segment *> &segs,
  int (*cmp)(const segment *, const segment *), int low, int high)
{
  segment *pivot = segs[high];

  int i = low-1;

  for (int j = low; j < high; j++) {
    if (cmp(segs[j], pivot) <= 0) {
      i++;
      std::swap(segs[i], segs[j]);
    }
  }
  i++;
  std::swap(segs[i], segs[high]);
  return i;
}

static
void _quicksort_segments(std::vector<segment *> &segs,
  int (*cmp)(const segment *, const segment *), int start, int end)
{
  if (start >= end) return;

  int p = _quicksort_partition(segs, cmp, start, end);

  _quicksort_segments(segs, cmp, start, p-1);
  _quicksort_segments(segs, cmp, p+1, end);
}

static
void _report_intersection(point *p,
    const std::vector<segment *> Lp,
    const std::vector<segment *> Cp,
    const std::vector<segment *> Up)
{
    printf("    INTERSECTION: id: P%lu - (%f %f)\n", p->id, p->x, p->y);
    printf("                  Ending: ");
    for (segment *s : Lp)
        printf("S%lu ", s->id);
    puts("");
    printf("                  Passing: ");
    for (segment *s : Cp)
        printf("S%lu ", s->id);
    puts("");
    printf("                  Starting: ");
    for (segment *s : Up)
        printf("S%lu ", s->id);
    puts("");
}

static
void _find_new_event(event *&Q, point *p, status *S0, status *S1,
    std::vector<point *> &points)
{
    // Do the segments intersect?
    segment *s0 = S0->s;
    segment *s1 = S1->s;

    point *a = s0->p;
    point *b = s0->q;
    point *c = s1->p;
    point *d = s1->q;

    E_Float denom = a->x * (d->y - c->y) +
                   b->x * (c->y - d->y) +
                   d->x * (b->y - a->y) +
                   c->x * (a->y - b->y);

    // Parallel / Collinear: skip for now
    if (get_sign(denom) == 0) return;

    E_Float num = a->x * (d->y - c->y) +
                 c->x * (a->y - d->y) +
                 d->x * (c->y - a->y);

    E_Float s = num / denom;

    if (s < -TOL || s > 1.0+TOL) return;

    num = - ( a->x * (c->y - b->y) +
              b->x * (a->y - c->y) +
              c->x * (b->y - a->y) );

    E_Float t = num / denom;

    if (t < -TOL || t > 1.0+TOL) return;

    E_Float x = a->x + s * (b->x - a->x);
    E_Float y = a->y + s * (b->y - a->y);

    // Where does the intersection lie wrt to sweep point
    int cmp = cmp_xyz(p->x, p->y, x, y);

    if (cmp == 1) return;

    // Intersection is to the right of sweep point

    // Does a point with the same coordinate already exist in Q?
    event *I = event_locate(Q, x, y);

    if (I == NULL) {
        // No: insert new event and set its sit to s0
        point *P = new point(x, y, points.size(), s0->id);
        points.push_back(P);
        I = event_insert(Q, P);
    }

    // Set intersection info of S1 to I
    I->s = s0;
    S0->p = I->p;
}

static
void _handle_event(event *E, event *&Q, status *&T,
    const std::vector<segment *> &segs, size_t &seg_start,
    E_Float &xs, E_Float &ys,
    std::vector<point *> &points,
    std::unordered_map<point *, std::vector<segment *>> &Xsegs)
{
    // Handle passing and ending segments
    segment *s = E->s;
    point *p = E->p;

    std::vector<segment *> Up;
    
    while (seg_start < segs.size() && segs[seg_start]->p == p) {
        Up.push_back(segs[seg_start]);
        seg_start++;
    }

    status *sit = NULL;

    if (s == NULL) {
        segment sdummy(p, p, -1, -1);
        sit = status_locate(T, &sdummy, xs, ys);
        if (sit)
            s = sit->s;
    }

    std::vector<segment *> Lp;
    std::vector<segment *> Cp;
    
    // p is an intersection between Up, Lp and Cp

    segment *s_succ = NULL;
    segment *s_pred = NULL;

    if (s != NULL) {
        sit = status_lookup(T, s, xs, ys);
        
        // Walk up until you reach sit whose xit is either null or not p
        while (sit->p == p || (sit->c != NULL && (sit->c == status_succ(T, sit->s, xs, ys)->s)))
            sit = status_succ(T, sit->s, xs, ys);

        status *sit_succ = status_succ(T, sit->s, xs, ys);
        s_succ = sit_succ->s;

        // Walk down
        do {
            s = sit->s;

            if (s->q == p) {
                // Ending segment
                Lp.push_back(s);

                // Segment is to be deleted, prev segment no longer is collinear with it
                status *sit_prev = status_pred(T, s, xs, ys);
                if (sit_prev->c == s) {
                    sit_prev->c = NULL;
                    sit_prev->p = s->q;
                } 

                sit = sit_prev;
            } else {
                // If s is not collinear with its successor, change its intersection info to null
                // Will be done in later section
                Cp.push_back(s);
                sit = status_pred(T, s, xs, ys);
            }

        } while (sit->p == p || (sit->c != NULL && (sit->c == status_succ(T, sit->s, xs, ys)->s)));

        s_pred = sit->s;
    }

    // Update intersecting segments
    for (segment *s : Lp) Xsegs[p].push_back(s);
    for (segment *s : Cp) Xsegs[p].push_back(s);
    for (segment *s : Up) Xsegs[p].push_back(s);

    //if (Lp.size() + Cp.size() + Up.size() > 1)
    //    _report_intersection(p, Lp, Cp, Up);

    for (segment *s : Lp)
        T = status_delete(T, s, xs, ys);

    // Before advancing the sweep line, cache col info
    std::vector<segment *> col_info(Cp.size());

    for (size_t i = 0; i < Cp.size(); i++) {
        sit = status_lookup(T, Cp[i], xs, ys);
        col_info[i] = sit->c;
        T = status_delete(T, Cp[i], xs, ys);
    }

    assert(col_info.size() == Cp.size());

    xs = p->x;
    ys = p->y;

    for (size_t i = 0; i < col_info.size(); i++) {
        sit = status_insert(T, Cp[i], xs, ys);
        sit->c = col_info[i];
    }

    // get_succ(k): returns status<k'> where k' > k, null otherwise
    // locate_succ(k): returns status<k'> where k' >= 0, null otherwise

    for (size_t i = 0; i < Up.size(); i++) {
        segment *next_seg = Up[i]; 

        status_insert(T, next_seg, xs, ys); 
        event *xit = event_lookup(Q, next_seg->q);
        xit->s = next_seg;
        
        status *s_sit = status_succ(T, next_seg, xs, ys);
        status *p_sit = status_pred(T, next_seg, xs, ys); 

        if (s_succ == NULL) {
            s_succ = status_succ(T, next_seg, xs, ys)->s;
            s_pred = status_pred(T, next_seg, xs, ys)->s;
        }

        // succ segment should never be collinear
        assert(!segments_are_colli(s_sit->s, next_seg));

        // Attach p_seg to next_seg if they overlap

        s = p_sit->s;

        if (segments_are_colli(next_seg, s))
            p_sit->c = next_seg;
    }

    // Find new events
    if (s_succ != NULL) {
        status *sit_succ = status_lookup(T, s_succ, xs, ys);
        status *sit_pred = status_lookup(T, s_pred, xs, ys);

        status *sit_first = status_succ(T, sit_pred->s, xs, ys);
        _find_new_event(Q, p, sit_pred, sit_first, points);
        
        status *sit_last = status_pred(T, sit_succ->s, xs, ys);
        if (sit_last != sit_pred)
            _find_new_event(Q, p, sit_last, sit_succ, points);
    }
}

void _sweep(std::vector<segment *> &segs, std::vector<point *> &points,
    std::unordered_map<point *, std::vector<segment *>> & Xsegs)
{
    // Make event queue
    event *Q = NULL;
    
    for (size_t i = 0; i < segs.size(); i++) {
        segment *s = segs[i];
        assert(segment_is_lexico(s));
        event_insert(Q, s->p);
        event_insert(Q, s->q);
    }

    _quicksort_segments(segs, cmp_segments_lexico, 0, segs.size()-1);

    for (size_t i = 0; i < segs.size()-1; i++)
      assert(cmp_segments_lexico(segs[i], segs[i+1]) <= 0);

    //std::reverse(segs.begin(), segs.end());

    // Upper bound for the input coordinates
    E_Float BIG = 1;

    for (size_t i = 0; i < points.size(); i++) {
        BIG = std::max(fabs(points[i]->x), BIG);
        BIG = std::max(fabs(points[i]->y), BIG);
    }

    BIG *= 2;

    // Initialize status
    status *T = NULL;
    E_Float xs, ys;
    xs = ys = -BIG;

    // Add sentinels
    point *pinf0 = new point(-BIG, -BIG, -1, -1);
    point *pinf1 = new point( BIG, -BIG, -2, -2);
    point *pinf2 = new point(-BIG,  BIG, -3, -3);
    point *pinf3 = new point( BIG,  BIG, -4, -4);
    segment *sinf0 = new segment(pinf0, pinf1, segs.size(), segs.size());
    segment *sinf1 = new segment(pinf2, pinf3, segs.size()+1, segs.size()+1);

    status_insert(T, sinf0, xs, ys);
    status_insert(T, sinf1, xs, ys);

    // Go !
    size_t seg_start = 0;

    while (Q != NULL) {
        event *E = event_min(Q);
        //printf("Queue:\n");
        //event_print(Q);
        
        _handle_event(E, Q, T, segs, seg_start, xs, ys, points, Xsegs);
        
        Q = event_delete(Q, E->p);
        
        //printf("Status:\n");
        //status_print(T);
        //puts("");
    }

    assert(seg_start == segs.size());

    for (const auto P : Xsegs) {
      point *p = P.first;
      const auto &segs = P.second;
      printf("P%lu (%f %f): ", p->id, p->x, p->y);
      for (segment *s : segs)
        printf("S%lu ", s->id);
      puts("");
    }

    // Clean up
    delete pinf0;
    delete pinf1;
    delete pinf2;
    delete pinf3;
    delete sinf0;
    delete sinf1;
}
