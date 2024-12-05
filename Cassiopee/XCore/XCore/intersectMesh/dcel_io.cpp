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
#include "dcel.h"

void Dcel::write_face(const char *fname, const Face *f) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    E_Int np = 1;
    Hedge *h = f->rep;
    Hedge *w = h->next;
    while (w != h) {
        np++;
        w = w->next;
    }
    
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", np);

    h = f->rep;
    Vertex *v = h->orig;
    fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    w = h->next;
    while (w != h) {
        v = w->orig;
        fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
        w = w->next;
    }
    fclose(fh);
}

void Dcel::write_faces(const char *fname, const std::vector<Face *> &faces,
    E_Float scale) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    E_Int NP = 0;
    for (const Face *f : faces) {
        E_Int np = 1;
        Hedge *h = f->rep;
        Hedge *w = h->next;
        while (w != h) {
            np++;
            w = w->next;
        }
        NP += np;
    }

    fprintf(fh, "POINTS\n");
    fprintf(fh, "%d\n", NP);

    for (const Face *f : faces) {
        Hedge *h = f->rep;
        Vertex *v = h->orig;
        fprintf(fh, "%f %f %f\n", scale*v->x, scale*v->y, scale*v->z);
        Hedge *w = h->next;
        while (w != h) {
            v = w->orig;
            fprintf(fh, "%f %f %f\n", scale*v->x, scale*v->y, scale*v->z);
            w = w->next;
        }
    }

    fclose(fh);
}

void Dcel::write_hedge(const char *fname, const Hedge *h) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "2\n");
    Vertex *p = h->orig;
    Vertex *q = h->twin->orig;
    fprintf(fh, "%f %f %f\n", p->x, p->y, p->z);
    fprintf(fh, "%f %f %f\n", q->x, q->y, q->z);
    fprintf(fh, "EDGES\n");
    fprintf(fh, "1\n");
    fprintf(fh, "0 1\n");
    fclose(fh);
}

void Dcel::write_vertex(const char *fname, const Vertex *v) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "1\n");
    fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    fclose(fh);
}

void Dcel::write_point(const char *fname,
    const std::vector<Dcel::Vertex *> &I) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);
    fprintf(fh, "POINTS\n");
    fprintf(fh, "%zu\n", I.size());
    for (auto &v : I) fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    fclose(fh);
}

void Dcel::write_cycles_of_type(const char *fname, int type) const
{
    auto indices = extract_indices_of_type(type);
    auto cids = extract_cycles_of_indices(indices);
    write_ngon(fname, cids);
}

void Dcel::write_degen_cycles(const char *fname) const
{
    write_cycles_of_type(fname, Cycle::DEGEN);
}

void Dcel::write_inner_cycles(const char *fname) const
{
    write_cycles_of_type(fname, Cycle::INNER);
}

void Dcel::write_hole_cycles(const char *fname) const
{
    write_cycles_of_type(fname, Cycle::HOLE);
}

void Dcel::write_outer_cycles(const char *fname) const
{
    write_cycles_of_type(fname, Cycle::OUTER);
}

void Dcel::write_ngon(const char *fname, const std::vector<Cycle *> &cycles) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    E_Int np = 0;
    E_Int ne = 0;
    E_Int nf = (E_Int)cycles.size();

    std::map<Vertex *, E_Int> vmap;
    std::vector<Vertex *> new_pids;

    for (size_t i = 0; i < cycles.size(); i++) {
        Cycle *c = cycles[i];
        Hedge *h = c->rep;
        ne++;
        Vertex *p = h->orig;
        if (vmap.find(p) == vmap.end()) {
            vmap[p] = np++;
            new_pids.push_back(p);
        }
        Hedge *w = h->next;
        while (w != h) {
            p = w->orig;
            if (vmap.find(p) == vmap.end()) {
                vmap[p] = np++;
                new_pids.push_back(p);
            }
            ne++;
            w = w->next;
        }
    }

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (const auto &v : new_pids) {
        fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    }
    
    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", ne+1);
    E_Int sizeNGon = 0;
    fprintf(fh, SF_D_ " ", sizeNGon);
    for (E_Int i = 0; i < ne; i++) {
        sizeNGon += 2;
        fprintf(fh, SF_D_ " ", sizeNGon);
    }
    assert(sizeNGon == 2*ne);
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, SF_D_ "\n", sizeNGon);
    for (Cycle *c : cycles) {
        Hedge *h = c->rep;
        Vertex *p = h->orig;
        Vertex *q = h->twin->orig;
        fprintf(fh, SF_D_ " "  SF_D_ " ", vmap[p], vmap[q]);
        Hedge *w = h->next;
        while (w != h) {
            p = w->orig;
            q = w->twin->orig;
            fprintf(fh, SF_D_ " " SF_D_ " ", vmap[p], vmap[q]);
            w = w->next;
        }
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    E_Int sizeNFace = 0;
    fprintf(fh, SF_D_ " ", sizeNFace);
    for (Cycle *c : cycles) {
        Hedge *h = c->rep;
        sizeNFace += 1;
        Hedge *w = h->next;
        while (w != h) {
            sizeNFace += 1;
            w = w->next;
        }
        fprintf(fh, SF_D_ " ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, SF_D_ "\n", sizeNFace);
    for (E_Int i = 0; i < sizeNFace; i++)
        fprintf(fh, SF_D_ " ", i);

    fclose(fh);
}

void Dcel::write_ngon(const char *fname, const std::vector<E_Int> &fids) const
{
    std::vector<Face *> faces;
    for (E_Int fid : fids) faces.push_back(F[fid]);
    write_ngon(fname, faces);
}

void Dcel::write_ngon(const char *fname, const std::vector<Face *> &faces) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    E_Int np = 0;
    E_Int ne = 0;
    E_Int nf = (E_Int)faces.size();

    std::map<Vertex *, E_Int> vmap;
    std::vector<Vertex *> new_pids;

    for (size_t i = 0; i < faces.size(); i++) {
        Face *c = faces[i];
        Hedge *h = c->rep;
        ne++;
        Vertex *p = h->orig;
        if (vmap.find(p) == vmap.end()) {
            vmap[p] = np++;
            new_pids.push_back(p);
        }
        Hedge *w = h->next;
        while (w != h) {
            p = w->orig;
            if (vmap.find(p) == vmap.end()) {
                vmap[p] = np++;
                new_pids.push_back(p);
            }
            ne++;
            w = w->next;
        }
    }

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (const auto &v : new_pids) {
        fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    }

    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", ne+1);
    E_Int sizeNGon = 0;
    fprintf(fh, SF_D_ " ", sizeNGon);
    for (E_Int i = 0; i < ne; i++) {
        sizeNGon += 2;
        fprintf(fh, SF_D_ " ", sizeNGon);
    }
    assert(sizeNGon == 2*ne);
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, SF_D_ "\n", sizeNGon);
    for (Face *c : faces) {
        Hedge *h = c->rep;
        Vertex *p = h->orig;
        Vertex *q = h->twin->orig;
        fprintf(fh, SF_D_ " "  SF_D_ " ", vmap[p], vmap[q]);
        Hedge *w = h->next;
        while (w != h) {
            p = w->orig;
            q = w->twin->orig;
            fprintf(fh, SF_D_ " " SF_D_ " ", vmap[p], vmap[q]);
            w = w->next;
        }
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    E_Int sizeNFace = 0;
    fprintf(fh, SF_D_ " ", sizeNFace);
    for (Face *c : faces) {
        Hedge *h = c->rep;
        sizeNFace += 1;
        Hedge *w = h->next;
        while (w != h) {
            sizeNFace += 1;
            w = w->next;
        }
        fprintf(fh, SF_D_ " ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, SF_D_ "\n", sizeNFace);
    for (E_Int i = 0; i < sizeNFace; i++)
        fprintf(fh, SF_D_ " ", i);

    fclose(fh);
}