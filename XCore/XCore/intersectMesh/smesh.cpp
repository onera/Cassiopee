#include "proto.h"
#include <queue>

smesh::smesh(Mesh *M, E_Int *patch, E_Int nf)
{
  F.resize(nf);
  for (E_Int i = 0; i < nf; i++) {
    E_Int face = patch[i];
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, M);
    assert(np == 4 || np == 3);
    auto &cn = F[i];
    cn.resize(np);
    for (E_Int j = 0; j < np; j++) cn[j] = pn[j];
  }

  E_Int np = 0;

  for (size_t i = 0; i < F.size(); i++) {
    auto &cn = F[i];
    for (size_t j = 0; j < cn.size(); j++) {
      E_Int gp = cn[j];
      auto it = g2l_p.find(gp);
      if (it == g2l_p.end()) {
        // New point
        g2l_p[gp] = np;
        np++;
      }
    }
  }

  printf("Smesh points: %d\n", np);

  // Make X Y
  X.resize(np);
  Y.resize(np);

  for (size_t i = 0; i < F.size(); i++) {
    auto &cn = F[i];
    for (size_t j = 0; j < cn.size(); j++) {
      E_Int gp = cn[j];
      E_Int lp = g2l_p[gp];
      X[lp] = M->x[gp];
      Y[lp] = M->y[gp];
    }
  }

  // Replace with local ids
  for (auto &cn : F) {
    for (auto &p : cn) {
      p = g2l_p[p];
      assert(p >= 0 && p < np);
    }
  }

  make_edges();
}

o_edge::o_edge(E_Int P, E_Int Q)
: p(P), q(Q)
{}

struct o_edge_cmp {
  bool operator()(const o_edge &e, const o_edge &f) const
  {
    E_Int e_p = std::min(e.p, e.q);
    E_Int e_q = std::max(e.p, e.q);
    E_Int f_p = std::min(f.p, f.q);
    E_Int f_q = std::max(f.p, f.q);
    return (e_p < f_p) || (e_p == f_p && e_q < f_q);
  }
};

void smesh::make_edges()
{
  // Order the faces counter-clockwise
  for (size_t i = 0; i < F.size(); i++) {
    auto &face = F[i];
    E_Int a = face[0];
    E_Int b = face[1];
    E_Int c = face[2];
    E_Float ax = X[a];
    E_Float ay = Y[a];
    E_Float bx = X[b];
    E_Float by = Y[b];
    E_Float cx = X[c];
    E_Float cy = Y[c];
    E_Float det = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
    int s = sign(det);
    assert(s);
    if (s < 0) {
      std::reverse(face.begin(), face.end());
    }
  }

  // Make the edges
  F2E.resize(F.size());
  std::map<o_edge, int, o_edge_cmp> edges;

  for (size_t i = 0; i < F.size(); i++) {
    auto &face = F[i];
    for (size_t j = 0; j < face.size(); j++) {
      E_Int p = face[j];
      E_Int q = face[(j+1)%face.size()];
      o_edge e(p, q);
      auto it = edges.find(e);
      if (it == edges.end()) {
        // New edge
        F2E[i].push_back(E.size());
        edges[e] = E.size();
        E.push_back(e);
      } else {
        F2E[i].push_back(it->second);
      }
    }
  }

  // Euler's formula should always be true
  assert(F.size()+1 + X.size() == E.size() + 2);

  // Make edge_to_face
  E2F.resize(E.size(), {-1, -1});
  for (size_t i = 0; i < F2E.size(); i++) {
    auto &face = F2E[i];
    for (size_t j = 0; j < face.size(); j++) {
      E_Int e = face[j];
      if (E2F[e][0] == -1) E2F[e][0] = i;
      else E2F[e][1] = i;
    }
  }

  // Make faces neighborhood
  std::vector<std::vector<E_Int>> F2F(F.size());
  for (size_t i = 0; i < F2E.size(); i++) {
    auto &face = F2E[i];
    auto &neis = F2F[i];
    for (size_t j = 0; j < face.size(); j++) {
      E_Int e = face[j];
      if (E2F[e][0] == (E_Int)i) neis.push_back(E2F[e][1]);
      else if (E2F[e][1] == (E_Int)i) neis.push_back(E2F[e][0]);
      else assert(0);
    }
  }

  // Traverse the facelist breadth-first and adjust edges accordingly
  std::vector<int> visited(F.size(), 0);
  std::queue<E_Int> Q;
  Q.push(0);
  visited[0] = 1;

  while (!Q.empty()) {
    E_Int f = Q.front();
    Q.pop();

    assert(f != -1);

    visited[f] = 1;

    auto &neis = F2F[f];
    auto &edgs = F2E[f];
    auto &face = F[f];

    for (size_t j = 0; j < face.size(); j++) {
      E_Int nei = neis[j];

      E_Int p = face[j];
      E_Int q = face[(j+1)%face.size()];

      E_Int e = edgs[j];

      if (nei == -1) {
        assert(E[e].p == p);
        assert(E[e].q == q);
        continue;
      }

      if (visited[nei]) {
        assert(E2F[e][0] == nei);
        assert(E2F[e][1] == f);
        assert(E[e].p == q);
        assert(E[e].q == p);
        continue;
      }

      if (E[e].p != p) {
        assert(visited[nei] == 0);
        assert(E[e].q == p);
        assert(E[e].p == q);
        std::swap(E[e].p, E[e].q);
        E2F[e][0] = f;
        E2F[e][1] = nei;
        Q.push(nei);
      }
    }
  }

  // Check everything is well set
  for (size_t i = 0; i < F.size(); i++) {
    auto &face = F[i];
    for (size_t j = 0; j < face.size(); j++) {
      E_Int e = F2E[i][j];
      E_Int p = face[j];
      E_Int q = face[(j+1)%face.size()];

      if (E[e].p == p) {
        assert(E[e].q == q);
        assert(E2F[e][0] == (E_Int)i);
      } else if (E[e].q == p) {
        assert(E[e].p == q);
        assert(E2F[e][1] == (E_Int)i);
      } else {
        assert(0);
      }
    }
  }

  puts("EDGES OK.");
}