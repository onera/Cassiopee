#include "xcore.h"
#include "common/Karray.h"
#include <set>
#include <unordered_set>

#define TOL 1e-12

struct Vtx {
    E_Float x, y, z;
    E_Int i0;
    mutable E_Int i1;
    E_Int id;

    Vtx(E_Float x_, E_Float y_, E_Float z_)
    : x(x_), y(y_), z(z_), i0(-1), i1(-1), id(-1)
    {}

    bool operator<(const Vtx &v) const
    {
        if (x < v.x) return true;
        if (x == v.x && y < v.y) return true;
        if (x == v.x && y == v.y && z < v.z) return true;
        return false;
    }
};

struct FACE {
    std::vector<E_Int> F;
    E_Int i0;
    mutable E_Int i1;
    E_Int id;

    FACE(E_Int *pn, E_Int np)
    {
        F.reserve(np);
        for (E_Int i = 0; i < np; i++) {
            F.push_back(pn[i]);
        }
        std::sort(F.begin(), F.end());
        i0 = i1 = id = -1;
    }

    bool operator==(const FACE &face) const
    {
        if (F.size() == face.F.size()) {
            for (size_t i = 0; i < F.size(); i++) {
                if (F[i] != face.F[i]) {
                    return false;
                }
            }

            return true;
        }

        return false;
    }
};

struct FACEHASH {
    uint64_t hash(uint64_t val, uint64_t seed) const
    {
        uint64_t HASH = seed;
        HASH += val;
        HASH += HASH << 10;
        HASH ^= HASH >> 6;
        return HASH;
    }

    uint32_t operator()(const FACE &face) const
    {
        uint64_t res = 0;
        for (size_t i = 0; i < face.F.size(); i++)
            res = hash(face.F[i], res);
      
        res += res << 3;
        res ^= res >> 11;
        res += res << 15;

        return res;
    }
};

PyObject *K_XCORE::IntersectMesh_Merge(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE;
    if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVE)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray, sarray;
    E_Int ret = Karray_parse_ngon(MASTER, marray);
    if (ret != 0) {
        RAISE("First input mesh is not an NGon.");
        return NULL;
    }

    ret = Karray_parse_ngon(SLAVE, sarray);
    if (ret != 0) {
        RAISE("Second input mesh is not an NGon.");
        Karray_free_ngon(marray);
        return NULL;
    }

    // Points
    std::set<Vtx> points;

    // Master points
    E_Float *X = marray.x;
    E_Float *Y = marray.y;
    E_Float *Z = marray.z;

    printf("M points: %d\n", marray.npts);
    printf("S points: %d\n", sarray.npts);

    E_Int NP = 0;

    for (E_Int i = 0; i < marray.npts; i++) {
        Vtx xyz(X[i], Y[i], Z[i]);
        auto it = points.find(xyz);
        if (it != points.end()) abort();
        xyz.i0 = i;
        xyz.id = NP;
        points.insert(xyz);
        NP++;
    }

    assert((size_t)NP == points.size());

    // Slave points
    X = sarray.x;
    Y = sarray.y;
    Z = sarray.z;

    for (E_Int i = 0; i < sarray.npts; i++) {
        Vtx xyz(X[i], Y[i], Z[i]);
        auto it = points.find(xyz);
        if (it == points.end()) {
            xyz.i1 = i;
            xyz.id = NP;
            points.insert(xyz);
            NP++;
        } else {
            assert(it->i0 != -1);
            assert(it->id != -1);
            it->i1 = i;
        }
    }

    printf("Duplicate points: %lu\n", (size_t)(marray.npts + sarray.npts) - points.size());

    // Change the points ids within marray

    X = marray.x;
    Y = marray.y;
    Z = marray.z;

    for (E_Int fid = 0; fid < marray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = marray.get_face(fid, np);
        for (E_Int j = 0; j < np; j++) {
            E_Int pid = pn[j] - 1;
            Vtx xyz(X[pid], Y[pid], Z[pid]);
            auto it = points.find(xyz);
            if (it == points.end()) abort();
            pn[j] = it->id;
        }
    }

    // Change the points within sarray

    X = sarray.x;
    Y = sarray.y;
    Z = sarray.z;

    for (E_Int fid = 0; fid < sarray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = sarray.get_face(fid, np);
        for (E_Int j = 0; j < np; j++) {
            E_Int pid = pn[j] - 1;
            Vtx xyz(X[pid], Y[pid], Z[pid]);
            auto it = points.find(xyz);
            if (it == points.end()) abort();
            pn[j] = it->id;
        }
    }

    // Hash the faces
    std::unordered_set<FACE, FACEHASH> faces;
    
    // Hash the M faces

    E_Int NF = 0;

    for (E_Int fid = 0; fid < marray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = marray.get_face(fid, np);
        FACE face(pn, np);
        auto it = faces.find(face);
        if (it != faces.end()) abort();
        face.i0 = fid;
        face.id = NF;
        faces.insert(face);
        NF++; 
    }

    // Hash the S faces

    for (E_Int fid = 0; fid < sarray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = sarray.get_face(fid, np);
        FACE face(pn, np);
        auto it = faces.find(face);
        if (it == faces.end()) {
            face.i1 = fid;
            face.id = NF;
            faces.insert(face);
            NF++;
        } else {
            assert(it->i0 != -1);
            assert(it->id != -1);
            it->i1 = fid;
        }
    }

    printf("Duplicate faces: %lu\n", marray.nfaces() + sarray.nfaces() - faces.size());

    // Replace face ids in marray

    for (E_Int cid = 0; cid < marray.ncells(); cid++) {
        E_Int nf = -1;
        E_Int *pf = marray.get_cell(cid, nf);
        for (E_Int i = 0; i < nf; i++) {
            E_Int fid = pf[i]-1;
            E_Int np = -1;
            E_Int *pn = marray.get_face(fid, np);
            FACE face(pn, np);
            auto it = faces.find(face);
            if (it == faces.end()) abort();
            pf[i] = it->id;
        }
    }


    // Replace face ids in sarray

    for (E_Int cid = 0; cid < sarray.ncells(); cid++) {
        E_Int nf = -1;
        E_Int *pf = sarray.get_cell(cid, nf);
        for (E_Int i = 0; i < nf; i++) {
            E_Int fid = pf[i]-1;
            E_Int np = -1;
            E_Int *pn = sarray.get_face(fid, np);
            FACE face(pn, np);
            auto it = faces.find(face);
            if (it == faces.end()) abort();
            pf[i] = it->id;
        }
    }

    
    E_Int sizeNGon = 0, sizeNFace = 0;

    for (E_Int cid = 0; cid < marray.ncells(); cid++) {
        E_Int nf = -1;
        marray.get_cell(cid, nf);
        sizeNFace += nf;
    }

    for (E_Int cid = 0; cid < sarray.ncells(); cid++) {
        E_Int nf = -1;
        sarray.get_cell(cid, nf);
        sizeNFace += nf;
    }

    printf("sizeNFace: %d\n", sizeNFace);

    for (E_Int fid = 0; fid < marray.nfaces(); fid++) {
        E_Int np = -1;
        marray.get_face(fid, np);
        sizeNGon += np;
    }

    for (E_Int fid = 0; fid < sarray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = sarray.get_face(fid, np);
        // Don't add it if it's a dup
        FACE face(pn, np);
        auto it = faces.find(face);
        if (it->i0 == -1)
            sizeNGon += np;
    }

    printf("sizeNGon: %d\n", sizeNGon);

    const char *var_string = "CoordinateX,CoordinateY,CoordinateZ";

    E_Int NC = marray.ncells() + sarray.ncells();

    printf("NC: %d\n", NC);
    printf("NF: %d\n", NF);
    printf("NP: %d\n", NP);

    PyObject *array = K_ARRAY::buildArray3(3, var_string, NP, NC, NF,
        "NGON", sizeNGon, sizeNFace, 3, false, 3);

    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(array, f, cn);

    E_Float *px = f->begin(1);
    E_Float *py = f->begin(2);
    E_Float *pz = f->begin(3);

    for (const auto &point : points) {
        if (point.id >= NP) abort();
        if (point.id < 0) abort();
        px[point.id] = point.x;
        py[point.id] = point.y;
        pz[point.id] = point.z;
    }

    E_Int *indPG = cn->getIndPG();
    E_Int *ngon = cn->getNGon();
    E_Int *indPH = cn->getIndPH();
    E_Int *nface = cn->getNFace();


    indPG[0] = indPH[0] = 0;

    for (E_Int fid = 0; fid < marray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = marray.get_face(fid, np);
        FACE face(pn, np);
        auto it = faces.find(face);
        indPG[it->id+1] = np;
    }

    for (E_Int fid = 0; fid < sarray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = sarray.get_face(fid, np);
        FACE face(pn, np);
        auto it = faces.find(face);
        indPG[it->id+1] = np;
    }

    for (E_Int i = 0; i < NF; i++) indPG[i+1] += indPG[i];

    printf("indPG[NF]: %d\n", indPG[NF]);
    fflush(stdout);

    if (indPG[NF] != sizeNGon) abort();

    puts("ok indPG");
    fflush(stdout);

    for (E_Int cid = 0; cid < marray.ncells(); cid++) {
        E_Int nf = -1;
        marray.get_cell(cid, nf);
        indPH[cid+1] = nf;
    }

    for (E_Int cid = 0; cid < sarray.ncells(); cid++) {
        E_Int nf = -1;
        sarray.get_cell(cid, nf);
        indPH[cid + marray.ncells() + 1] = nf;
    }

    for (E_Int i = 0; i < NC; i++) indPH[i+1] += indPH[i];

    printf("indPH[NC]: %d\n", indPH[NC]);
    fflush(stdout);

    if (indPH[NC] != sizeNFace) abort();

    puts("ok indPH");
    fflush(stdout);

    for (E_Int fid = 0; fid < marray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = marray.get_face(fid, np);
        FACE face(pn, np);
        auto it = faces.find(face);
        E_Int *ptr = &ngon[indPG[it->id]];
        for (E_Int i = 0; i < np; i++)
            *ptr++ = pn[i] + 1;
    }

    for (E_Int fid = 0; fid < sarray.nfaces(); fid++) {
        E_Int np = -1;
        E_Int *pn = sarray.get_face(fid, np);
        FACE face(pn, np);
        auto it = faces.find(face);
        E_Int *ptr = &ngon[indPG[it->id]];
        for (E_Int i = 0; i < np; i++)
            *ptr++ = pn[i] + 1;
    }

    for (E_Int cid = 0; cid < marray.ncells(); cid++) {
        E_Int nf = -1;
        E_Int *pf = marray.get_cell(cid, nf);
        E_Int *ptr = &nface[indPH[cid]];
        for (E_Int i = 0; i < nf; i++)
            *ptr++ = pf[i] + 1;
    }

    for (E_Int cid = 0; cid < sarray.ncells(); cid++) {
        E_Int nf = -1;
        E_Int *pf = sarray.get_cell(cid, nf);
        E_Int *ptr = &nface[indPH[cid + marray.ncells()]];
        for (E_Int i = 0; i < nf; i++)
            *ptr++ = pf[i] + 1;
    }

    delete f;
    delete cn;

    return array;
}
