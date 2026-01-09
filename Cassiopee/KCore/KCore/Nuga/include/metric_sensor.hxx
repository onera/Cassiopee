/*    
    Copyright 2013-2026 ONERA.

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

//Authors : Sam Landier (sam.landier@onera.fr); Imad Hammani (imad.hammani@onera.fr)

#ifndef NUGA_METRIC_SENSOR_HXX
#define NUGA_METRIC_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/Metric.h"
//#include "Nuga/include/shell_smoother.hxx"
//#include "Nuga/include/V1_smoother.hxx"


namespace NUGA
{

template <typename mesh_t>
class metric_sensor : public sensor<mesh_t, DELAUNAY::VarMetric<DELAUNAY::Aniso3D>>
{
  public:
    bool _single_pass_done;
    
    //using cell_incr_t = DELAUNAY::Aniso3D;
    using sensor_input_t = DELAUNAY::VarMetric<DELAUNAY::Aniso3D>; // either vector<int> (ISO) or vector(incr_t<3>) (DIR)
    using parent_t = sensor<mesh_t, sensor_input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    metric_sensor(mesh_t& mesh) : parent_t(mesh, nullptr)
    {
        parent_t::_hmesh._sensor = this; /*in order to have an access to the field values from hmesh*/
        _canon_info.resize(mesh._ng.PHs.size());
        _start_nodes.resize(mesh._ng.PHs.size());
    }

    E_Int assign_data(const sensor_input_t& data) override;

    E_Int assign_data(const K_FLD::FloatArray& data);
    
    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override;

    void metric_fix(mesh_t& hmesh, output_t& adap_incr);

    void metric_fix_2(mesh_t& hmesh, output_t& adap_incr);

    bool stop() override { return _single_pass_done; }

    void Q4_adap_compute(E_Int PGi, output_t& adap_incr);

    void Hexa_adap_compute(E_Int PHi, output_t& adap_incr);

    bool is_FXY(E_Int PGi, output_t& adap_incr);

    bool is_HX12(E_Int PHi, output_t& adap_incr);

    bool is_HX18(E_Int PHi, output_t& adap_incr);

    bool is_HX27(E_Int PHi, output_t& adap_incr);

    std::vector<std::vector<E_Int>> _canon_info;
    std::vector<std::vector<E_Int>> _start_nodes;

    void compute_canon_info_bottom(E_Int PHi, E_Int n0);
    void compute_canon_info_top(E_Int PHi, E_Int n0);
    void compute_canon_info_left(E_Int PHi, E_Int n0);
    void compute_canon_info_right(E_Int PHi, E_Int n0);
    void compute_canon_info_front(E_Int PHi, E_Int n0);
    void compute_canon_info_back(E_Int PHi, E_Int n0);

    bool FIX_BOT_TOP(E_Int PHi, E_Int nei, output_t& adap_incr);
    bool FIX_LFT_RGT(E_Int PHi, E_Int nei, output_t& adap_incr);
    bool FIX_FRO_BCK(E_Int PHi, E_Int nei, output_t& adap_incr);

    void fix_faces(E_Int PHi, output_t& adap_incr, mesh_t& hmesh);

    E_Int prepare_bot_nei(E_Int PHi, mesh_t& hmesh, E_Int bot, E_Int i0, output_t& adap_incr);
    E_Int prepare_top_nei(E_Int PHi, mesh_t& hmesh, E_Int top, E_Int i0, output_t& adap_incr);
    E_Int prepare_lft_nei(E_Int PHi, mesh_t& hmesh, E_Int lft, E_Int i0, output_t& adap_incr);
    E_Int prepare_rgt_nei(E_Int PHi, mesh_t& hmesh, E_Int rgt, E_Int i0, output_t& adap_incr);
    E_Int prepare_fro_nei(E_Int PHi, mesh_t& hmesh, E_Int fro, E_Int i0, output_t& adap_incr);
    E_Int prepare_bck_nei(E_Int PHi, mesh_t& hmesh, E_Int bck, E_Int i0, output_t& adap_incr);

    void fix_alpha_beta(E_Int PHi, E_Int NEI, output_t& adap_incr);
    void fix_alpha_gamma(E_Int PHi, E_Int NEI, output_t& adap_incr);
    void fix_beta_gamma(E_Int PHi, E_Int NEI, output_t& adap_incr);

    bool enforce_consistancy(E_Int PHi, mesh_t& hmesh, output_t& adap_incr);
    void check_consistancy(E_Int PHi, mesh_t& hmesh, output_t& adap_incr);
};

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_bottom(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Bottom
    PGi = faces[0] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = n0;
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    nodes[0] = local[0];
    nodes[1] = local[1];
    nodes[2] = local[2];
    nodes[3] = local[3];
    _start_nodes[PHi][0] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    nodes[7] = local[2];
    nodes[4] = local[3];
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    _start_nodes[PHi][2] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    nodes[6] = local[2];
    nodes[5] = local[3];
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[1] == nodes[5]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[7]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    _start_nodes[PHi][1] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_top(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Top
    PGi = faces[1] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = n0;
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[4] = local[0];
    nodes[5] = local[1];
    nodes[6] = local[2];
    nodes[7] = local[3];
    _start_nodes[PHi][1] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[7]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[7]);
    assert(local[1] == nodes[4]);
    nodes[0] = local[2];
    nodes[3] = local[3];
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    _start_nodes[PHi][2] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[6]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[6]);
    assert(local[1] == nodes[5]);
    nodes[1] = local[2];
    nodes[2] = local[3];
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[1]);
    assert(local[2] == nodes[2]);
    assert(local[3] == nodes[3]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    _start_nodes[PHi][0] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_left(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Left
    PGi = faces[2] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = n0;
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[0] = local[0];
    nodes[3] = local[1];
    nodes[7] = local[2];
    nodes[4] = local[3];
    _start_nodes[PHi][2] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[3] == nodes[3]);
    nodes[1] = local[1];
    nodes[2] = local[2];
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[3] == nodes[7]);
    nodes[5] = local[1];
    nodes[6] = local[2];
    _start_nodes[PHi][1] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[5]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_right(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Right
    PGi = faces[3] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = n0;
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[1] = local[0];
    nodes[2] = local[1];
    nodes[6] = local[2];
    nodes[5] = local[3];
    _start_nodes[PHi][3] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[3] == nodes[5]);
    nodes[0] = local[1];
    nodes[4] = local[2];
    _start_nodes[PHi][4] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[3] == nodes[6]);
    nodes[3] = local[1];
    nodes[7] = local[2];
    _start_nodes[PHi][5] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[4]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    _start_nodes[PHi][2] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[1]);
    assert(local[2] == nodes[2]);
    assert(local[3] == nodes[3]);
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[1] == nodes[5]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[7]);
    _start_nodes[PHi][1] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_front(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Front
    PGi = faces[4] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = n0;
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[1] = local[0];
    nodes[0] = local[1];
    nodes[4] = local[2];
    nodes[5] = local[3];
    _start_nodes[PHi][4] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[1]);
    nodes[2] = local[2];
    nodes[3] = local[3];
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[4]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[4]);
    assert(local[1] == nodes[5]);
    nodes[6] = local[2];
    nodes[7] = local[3];
    _start_nodes[PHi][1] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[5]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[4]);
    _start_nodes[PHi][2] = i0;

    // Back
    PGi = faces[5] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[6]);
    _start_nodes[PHi][5] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::compute_canon_info_back(E_Int PHi, E_Int n0)
{
    auto& hmesh = parent_t::_hmesh;
    const auto& ng = hmesh._ng;
    const auto& F2E = hmesh._F2E;

    E_Int nodes[8];

    const E_Int *faces = ng.PHs.get_facets_ptr(PHi);

    bool reorient;
    E_Int i0;
    E_Int local[4];
    E_Int canon_data[12];
    E_Int PGi;

    _start_nodes[PHi].resize(6);

    // Back
    PGi = faces[5] - 1;
    const E_Int *pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = n0;
    canon_data[10] = (E_Int) reorient;
    canon_data[11] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == pN[i0]);
    nodes[2] = local[0];
    nodes[3] = local[1];
    nodes[7] = local[2];
    nodes[6] = local[3];
    _start_nodes[PHi][5] = i0;

    // Bottom
    PGi = faces[0] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[2]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[2]);
    assert(local[1] == nodes[3]);
    nodes[0] = local[2];
    nodes[1] = local[3];
    // local[2] is the usual i0 here
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[0] = (E_Int) reorient;
    canon_data[1] = i0;
    _start_nodes[PHi][0] = i0;

    // Top
    PGi = faces[1] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[6]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[6]);
    assert(local[1] == nodes[7]);
    nodes[4] = local[2];
    nodes[5] = local[3];
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, local[2]);
    canon_data[2] = (E_Int) reorient;
    canon_data[3] = i0;
    _start_nodes[PHi][1] = i0;

    // Right
    PGi = faces[3] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, false, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[2]);
    assert(local[2] == nodes[6]);
    assert(local[3] == nodes[5]);
    canon_data[6] = (E_Int) reorient;
    canon_data[7] = i0;
    _start_nodes[PHi][3] = i0;

    // Left
    PGi = faces[2] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[0]);
    canon_data[4] = (E_Int) reorient;
    canon_data[5] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[0]);
    assert(local[1] == nodes[3]);
    assert(local[2] == nodes[7]);
    assert(local[3] == nodes[4]);
    _start_nodes[PHi][2] = i0;

    // Front
    PGi = faces[4] - 1;
    pN = ng.PGs.get_facets_ptr(PGi);
    reorient = K_MESH::Quadrangle::need_a_reorient(PGi, PHi, true, F2E);
    i0 = K_CONNECT::IdTool::get_pos(pN, 4, nodes[1]);
    canon_data[8] = (E_Int) reorient;
    canon_data[9] = i0;
    K_MESH::Quadrangle::order_nodes(local, pN, reorient, i0);
    assert(local[0] == nodes[1]);
    assert(local[1] == nodes[0]);
    assert(local[2] == nodes[4]);
    assert(local[3] == nodes[5]);
    _start_nodes[PHi][4] = i0;

    _canon_info[PHi].resize(6);

    // deduce swap
    for (E_Int i = 0; i < 6; i++) {
        E_Int reorient = canon_data[2*i];
        E_Int i0 = canon_data[2*i+1];
        if (reorient) {
            if (i0%2 == 0) _canon_info[PHi][i] = 1;
            else _canon_info[PHi][i] = 0;
        } else {
          if (i0%2 == 0) _canon_info[PHi][i] = 0;
          else _canon_info[PHi][i] = 1;
        }
    }
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::assign_data(const K_FLD::FloatArray& arr) {
    // convert data to VarMetric<Aniso3D>

    E_Int npts = parent_t::_hmesh._crd.cols();

    assert(npts == arr.cols()); // input must be sized as number of pts (if decided to do differently,  need here to resize and fill missing field values)

    sensor_input_t data(parent_t::_hmesh._crd, -1, -1);
    DELAUNAY::Aniso3D m;

    // populate _data
    for (E_Int i = 0; i < npts; i++) {
        for (E_Int j = 0; j < 6; j++) {
            m[j] = arr(j,i);
        }
        if (!data.isValidMetric(m)) {
            fprintf(stderr, "BAD METRIC!\n");
            std::cout << "Point: " << i << std::endl;
            std::cout << m << std::endl;
            exit(1);
        }
        //assert(data.isValidMetric(m));
        data.setMetric(i,m);
    }

    std::cout << "assign_data ok\n";

    parent_t::_data = data;

  _single_pass_done = false;//reinit to enable

  return 0;

}


template <typename mesh_t>
bool metric_sensor<mesh_t>::is_FXY(E_Int PGi, output_t& adap_incr)
{
    if (adap_incr.face_adap_incr[PGi].n[0] > 0 && adap_incr.face_adap_incr[PGi].n[1] > 0)
        return true;
    return false;
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::is_HX12(E_Int PHi, output_t& adap_incr)
{
    E_Int count = 0;
    for (E_Int i = 0; i < 3; i++) {
        if (adap_incr.cell_adap_incr[PHi].n[i] > 0) count++;
    }

    return (count == 1);
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::is_HX18(E_Int PHi, output_t& adap_incr)
{
    E_Int count = 0;
    for (E_Int i = 0; i < 3; i++) {
        if (adap_incr.cell_adap_incr[PHi].n[i] > 0) count++;
    }

    return (count == 2);
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::is_HX27(E_Int PHi, output_t& adap_incr)
{
    E_Int count = 0;
    for (E_Int i = 0; i < 3; i++) {
        if (adap_incr.cell_adap_incr[PHi].n[i] > 0) count++;
    }

    return (count == 3);
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::prepare_bot_nei(E_Int PHi, mesh_t& hmesh, E_Int bot, E_Int i0, output_t& adap_incr)
{
    E_Int BOT = NEIGHBOR(PHi, hmesh._F2E, bot);
    if (BOT == IDX_NONE) return BOT;

    // make bot the top face of BOT
    E_Int *faces = hmesh._ng.PHs.get_facets_ptr(BOT);
    while (bot+1 != faces[1]) K_CONNECT::IdTool::right_shift<6>(faces, 1);

    // match the canon config of PHi
    K_MESH::Hexahedron::reorder_pgs_top(hmesh._ng, hmesh._F2E, BOT, i0);
    compute_canon_info_top(BOT, i0);
    Hexa_adap_compute(BOT, adap_incr);

    return BOT;
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::prepare_top_nei(E_Int PHi, mesh_t& hmesh, E_Int top, E_Int i0, output_t& adap_incr)
{
    E_Int TOP = NEIGHBOR(PHi, hmesh._F2E, top);
    if (TOP == IDX_NONE) return TOP;

    // make top the bot face of TOP
    E_Int *faces = hmesh._ng.PHs.get_facets_ptr(TOP);
    while (top+1 != faces[0]) K_CONNECT::IdTool::right_shift<6>(faces, 1);

    // match the canon config of PHi
    K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, TOP, i0);
    compute_canon_info_bottom(TOP, i0);
    Hexa_adap_compute(TOP, adap_incr);

    return TOP;
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::prepare_lft_nei(E_Int PHi, mesh_t& hmesh, E_Int lft, E_Int i0, output_t& adap_incr)
{
    E_Int LFT = NEIGHBOR(PHi, hmesh._F2E, lft);
    if (LFT == IDX_NONE) return LFT;

    // make lft the rgt face of LFT
    E_Int *faces = hmesh._ng.PHs.get_facets_ptr(LFT);
    while (lft+1 != faces[3]) K_CONNECT::IdTool::right_shift<6>(faces, 1);

    // match the canon config of PHi
    K_MESH::Hexahedron::reorder_pgs_right(hmesh._ng, hmesh._F2E, LFT, i0);
    compute_canon_info_right(LFT, i0);
    Hexa_adap_compute(LFT, adap_incr);

    return LFT;
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::prepare_rgt_nei(E_Int PHi, mesh_t& hmesh, E_Int rgt, E_Int i0, output_t& adap_incr)
{
    E_Int RGT = NEIGHBOR(PHi, hmesh._F2E, rgt);
    if (RGT == IDX_NONE) return RGT;

    // make rgt the lft face of RGT
    E_Int *faces = hmesh._ng.PHs.get_facets_ptr(RGT);
    while (rgt+1 != faces[2]) K_CONNECT::IdTool::right_shift<6>(faces, 1);

    // match the canon config of PHi
    K_MESH::Hexahedron::reorder_pgs_left(hmesh._ng, hmesh._F2E, RGT, i0);
    compute_canon_info_left(RGT, i0);
    Hexa_adap_compute(RGT, adap_incr);

    return RGT;
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::prepare_fro_nei(E_Int PHi, mesh_t& hmesh, E_Int fro, E_Int i0, output_t& adap_incr)
{
    E_Int FRO = NEIGHBOR(PHi, hmesh._F2E, fro);
    if (FRO == IDX_NONE) return FRO;

    // make fro the bck face of FRO
    E_Int *faces = hmesh._ng.PHs.get_facets_ptr(FRO);
    while (fro+1 != faces[5]) K_CONNECT::IdTool::right_shift<6>(faces, 1);

    // match the canon config of PHi
    K_MESH::Hexahedron::reorder_pgs_back(hmesh._ng, hmesh._F2E, FRO, i0);
    compute_canon_info_back(FRO, i0);
    Hexa_adap_compute(FRO, adap_incr);

    return FRO;
}

template <typename mesh_t>
E_Int metric_sensor<mesh_t>::prepare_bck_nei(E_Int PHi, mesh_t& hmesh, E_Int bck, E_Int i0, output_t& adap_incr)
{
    E_Int BCK = NEIGHBOR(PHi, hmesh._F2E, bck);
    if (BCK == IDX_NONE) return BCK;

    // make bck the fro face of BCK
    E_Int *faces = hmesh._ng.PHs.get_facets_ptr(BCK);
    while (bck+1 != faces[4]) K_CONNECT::IdTool::right_shift<6>(faces, 1);

    // match the canon config of PHi
    K_MESH::Hexahedron::reorder_pgs_front(hmesh._ng, hmesh._F2E, BCK, i0);
    compute_canon_info_front(BCK, i0);
    Hexa_adap_compute(BCK, adap_incr);

    return BCK;
}

template <typename mesh_t>
void metric_sensor<mesh_t>::fix_beta_gamma(E_Int PHi, E_Int NEI, output_t& adap_incr)
{
    E_Int incr[4] = {
        adap_incr.cell_adap_incr[PHi].n[1],
        adap_incr.cell_adap_incr[PHi].n[2],
        adap_incr.cell_adap_incr[NEI].n[1],
        adap_incr.cell_adap_incr[NEI].n[2]
    };

    E_Int INCR = 0;
    while (1) {
        if ((incr[0] == 0 && incr[1] == 0) || (incr[2] == 0 && incr[3] == 0)) break;

        else if (incr[0] == 0 && incr[2] == 0) break;

        else if (incr[1] == 0 && incr[3] == 0) break;

        else if (incr[0] == 0 && incr[1] == 0 && incr[2] == 0 && incr[3] == 0) break;

        // look for min
        E_Int MIN = IDX_NONE;
        //E_Int idx = -1;
        for (E_Int i = 0; i < 4; i++) 
        {
            if (incr[i] == 0) continue;
            if (incr[i] < MIN) 
            {
                MIN = incr[i];
                //idx = i;
            }
        }

        INCR += MIN;

        // decrement
        for (E_Int i = 0; i < 4; i++) 
        {
            incr[i] -= MIN;
            if (incr[i] < 0) incr[i] = 0;
        }
    }

    adap_incr.cell_adap_incr[PHi].n[1] = std::max(INCR, adap_incr.cell_adap_incr[PHi].n[1]);
    adap_incr.cell_adap_incr[PHi].n[2] = std::max(INCR, adap_incr.cell_adap_incr[PHi].n[2]);
    adap_incr.cell_adap_incr[NEI].n[1] = std::max(INCR, adap_incr.cell_adap_incr[NEI].n[1]);
    adap_incr.cell_adap_incr[NEI].n[2] = std::max(INCR, adap_incr.cell_adap_incr[NEI].n[2]);
}

template <typename mesh_t>
void metric_sensor<mesh_t>::fix_alpha_beta(E_Int PHi, E_Int NEI, output_t& adap_incr)
{
    E_Int incr[4] = {
        adap_incr.cell_adap_incr[PHi].n[0],
        adap_incr.cell_adap_incr[PHi].n[1],
        adap_incr.cell_adap_incr[NEI].n[0],
        adap_incr.cell_adap_incr[NEI].n[1]
    };

    E_Int INCR = 0;
    while (1) {
        if ((incr[0] == 0 && incr[1] == 0) || (incr[2] == 0 && incr[3] == 0)) break;

        else if (incr[0] == 0 && incr[2] == 0) break;

        else if (incr[1] == 0 && incr[3] == 0) break;

        else if (incr[0] == 0 && incr[1] == 0 && incr[2] == 0 && incr[3] == 0) break;

        // look for min
        E_Int MIN = IDX_NONE;
        E_Int idx = -1;
        for (E_Int i = 0; i < 4; i++) {
            if (incr[i] == 0) continue;
            if (incr[i] < MIN) {
                MIN = incr[i];
                idx = i;
            }
        }

        INCR += MIN;

        // decrement
        for (E_Int i = 0; i < 4; i++) {
            incr[i] -= MIN;
            if (incr[i] < 0) incr[i] = 0;
        }
    }

    adap_incr.cell_adap_incr[PHi].n[0] = std::max(INCR, adap_incr.cell_adap_incr[PHi].n[0]);
    adap_incr.cell_adap_incr[PHi].n[1] = std::max(INCR, adap_incr.cell_adap_incr[PHi].n[1]);
    adap_incr.cell_adap_incr[NEI].n[0] = std::max(INCR, adap_incr.cell_adap_incr[NEI].n[0]);
    adap_incr.cell_adap_incr[NEI].n[1] = std::max(INCR, adap_incr.cell_adap_incr[NEI].n[1]);

}

template <typename mesh_t>
void metric_sensor<mesh_t>::fix_alpha_gamma(E_Int PHi, E_Int NEI, output_t& adap_incr)
{
    E_Int incr[4] = {
        adap_incr.cell_adap_incr[PHi].n[0],
        adap_incr.cell_adap_incr[PHi].n[2],
        adap_incr.cell_adap_incr[NEI].n[0],
        adap_incr.cell_adap_incr[NEI].n[2]
    };

    E_Int INCR = 0;
    while (1) {
        if ((incr[0] == 0 && incr[1] == 0) || (incr[2] == 0 && incr[3] == 0)) break;

        else if (incr[0] == 0 && incr[2] == 0) break;

        else if (incr[1] == 0 && incr[3] == 0) break;

        else if (incr[0] == 0 && incr[1] == 0 && incr[2] == 0 && incr[3] == 0) break;

        // look for min
        E_Int MIN = IDX_NONE;
        E_Int idx = -1;
        for (E_Int i = 0; i < 4; i++) {
            if (incr[i] == 0) continue;
            if (incr[i] < MIN) {
                MIN = incr[i];
                idx = i;
            }
        }

        INCR += MIN;

        // decrement
        for (E_Int i = 0; i < 4; i++) {
            incr[i] -= MIN;
            if (incr[i] < 0) incr[i] = 0;
        }
    }

    adap_incr.cell_adap_incr[PHi].n[0] = std::max(INCR, adap_incr.cell_adap_incr[PHi].n[0]);
    adap_incr.cell_adap_incr[PHi].n[2] = std::max(INCR, adap_incr.cell_adap_incr[PHi].n[2]);
    adap_incr.cell_adap_incr[NEI].n[0] = std::max(INCR, adap_incr.cell_adap_incr[NEI].n[0]);
    adap_incr.cell_adap_incr[NEI].n[2] = std::max(INCR, adap_incr.cell_adap_incr[NEI].n[2]);
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::FIX_BOT_TOP(E_Int PHi, E_Int nei, output_t& adap_incr)
{
    if (nei == IDX_NONE) return false;

    // both cells are in coherent config
    E_Int ao = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int bo = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int co = adap_incr.cell_adap_incr[PHi].n[2];
    E_Int an = adap_incr.cell_adap_incr[nei].n[0];
    E_Int bn = adap_incr.cell_adap_incr[nei].n[1];
    E_Int cn = adap_incr.cell_adap_incr[nei].n[2];

    E_Int PHi_to_mod = -1;

    // fix alphas
    if (::abs(ao-an) <= 1) goto fix_beta;
    PHi_to_mod = ao > an ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[0] = std::max(ao, an) - 1;

    // fix betas
    fix_beta:
    if (::abs(bo-bn) <= 1) goto fix_gamma;
    PHi_to_mod = bo > bn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[1] = std::max(bo, bn) - 1;

    // fix gammas
    fix_gamma:
    if (::abs(co-cn) <= 1) goto FIX;
    PHi_to_mod = co > cn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[2] = std::max(co, cn) - 1;
    
    FIX:
    fix_alpha_beta(PHi, nei, adap_incr);

    bool has_changed = !(adap_incr.cell_adap_incr[PHi].n[0] == ao &&
             adap_incr.cell_adap_incr[PHi].n[1] == bo &&
             adap_incr.cell_adap_incr[PHi].n[2] == co);

    return has_changed;
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::FIX_LFT_RGT(E_Int PHi, E_Int nei, output_t& adap_incr)
{
    if (nei == IDX_NONE) return false;

    // both cells are in coherent config
    E_Int ao = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int bo = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int co = adap_incr.cell_adap_incr[PHi].n[2];
    E_Int an = adap_incr.cell_adap_incr[nei].n[0];
    E_Int bn = adap_incr.cell_adap_incr[nei].n[1];
    E_Int cn = adap_incr.cell_adap_incr[nei].n[2];

    E_Int PHi_to_mod = -1;

    // fix alphas
    if (::abs(ao-an) <= 1) goto fix_beta;
    PHi_to_mod = ao > an ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[0] = std::max(ao, an) - 1;

    // fix betas
    fix_beta:
    if (::abs(bo-bn) <= 1) goto fix_gamma;
    PHi_to_mod = bo > bn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[1] = std::max(bo, bn) - 1;

    // fix gammas
    fix_gamma:
    if (::abs(co-cn) <= 1) goto FIX;
    PHi_to_mod = co > cn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[2] = std::max(co, cn) - 1;

    FIX:
    fix_beta_gamma(PHi, nei, adap_incr);
    
    bool has_changed = !(adap_incr.cell_adap_incr[PHi].n[0] == ao &&
             adap_incr.cell_adap_incr[PHi].n[1] == bo &&
             adap_incr.cell_adap_incr[PHi].n[2] == co);

    return has_changed;
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::FIX_FRO_BCK(E_Int PHi, E_Int nei, output_t& adap_incr)
{
    if (nei == IDX_NONE) return false;

    // both cells are in coherent config
    E_Int ao = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int bo = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int co = adap_incr.cell_adap_incr[PHi].n[2];
    E_Int an = adap_incr.cell_adap_incr[nei].n[0];
    E_Int bn = adap_incr.cell_adap_incr[nei].n[1];
    E_Int cn = adap_incr.cell_adap_incr[nei].n[2];

    E_Int PHi_to_mod = -1;

    // fix alphas
    if (::abs(ao-an) <= 1) goto fix_beta;
    PHi_to_mod = ao > an ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[0] = std::max(ao, an) - 1;

    // fix betas
    fix_beta:
    if (::abs(bo-bn) <= 1) goto fix_gamma;
    PHi_to_mod = bo > bn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[1] = std::max(bo, bn) - 1;

    // fix gammas
    fix_gamma:
    if (::abs(co-cn) <= 1) goto FIX;
    PHi_to_mod = co > cn ? nei : PHi;
    adap_incr.cell_adap_incr[PHi_to_mod].n[2] = std::max(co, cn) - 1;
    
    FIX:
    fix_alpha_gamma(PHi, nei, adap_incr);
    
    bool has_changed = !(adap_incr.cell_adap_incr[PHi].n[0] == ao &&
             adap_incr.cell_adap_incr[PHi].n[1] == bo &&
             adap_incr.cell_adap_incr[PHi].n[2] == co);

    return has_changed;
}

template <typename mesh_t>
bool metric_sensor<mesh_t>::enforce_consistancy(E_Int PHi, mesh_t& hmesh, output_t& adap_incr)
{
    using face_incr_t = typename output_t::face_incr_t;

    E_Int alpha = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int beta = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int gamma = adap_incr.cell_adap_incr[PHi].n[2];

    const E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);
    const auto& swap = _canon_info[PHi];
    E_Int PGi = -1;

    face_incr_t f_incr[6];

    for (int i = 0; i < 6; i++) {
        PGi = pF[i] - 1;
        if (swap[i]) {
            f_incr[i].n[0] = adap_incr.face_adap_incr[PGi].n[1];
            f_incr[i].n[1] = adap_incr.face_adap_incr[PGi].n[0];
        } else {
            f_incr[i].n[0] = adap_incr.face_adap_incr[PGi].n[0];
            f_incr[i].n[1] = adap_incr.face_adap_incr[PGi].n[1];
        }
    }

    E_Int sillonX[4] = {f_incr[0].n[0], f_incr[4].n[0], f_incr[1].n[0], f_incr[5].n[0]};
    E_Int sillonY[4] = {f_incr[0].n[1], f_incr[2].n[0], f_incr[1].n[1], f_incr[3].n[0]};
    E_Int sillonZ[4] = {f_incr[2].n[1], f_incr[4].n[1], f_incr[3].n[1], f_incr[5].n[1]};

    adap_incr.cell_adap_incr[PHi].n[0] = std::max(adap_incr.cell_adap_incr[PHi].n[0], *std::max_element(sillonX, sillonX+4)-1);
    adap_incr.cell_adap_incr[PHi].n[1] = std::max(adap_incr.cell_adap_incr[PHi].n[1], *std::max_element(sillonY, sillonY+4)-1);
    adap_incr.cell_adap_incr[PHi].n[2] = std::max(adap_incr.cell_adap_incr[PHi].n[2], *std::max_element(sillonZ, sillonZ+4)-1);

    // return true if consistant
    return (adap_incr.cell_adap_incr[PHi].n[0] == alpha && 
            adap_incr.cell_adap_incr[PHi].n[1] == beta &&
            adap_incr.cell_adap_incr[PHi].n[2] == gamma);
}

template <typename mesh_t>
void metric_sensor<mesh_t>::metric_fix(mesh_t& hmesh, output_t& adap_incr)
{
    //using cell_incr_t = typename output_t::cell_incr_t;
    using face_incr_t = typename output_t::face_incr_t;

    std::stack<E_Int> stk;
    E_Int PGi = 0;

    for (; PGi < hmesh._ng.PGs.size(); PGi++) {
        if (adap_incr.face_adap_incr[PGi] != 0) stk.push(PGi);
    }

    E_Int own = -1, nei = -1;
    E_Int *pOwn = 0, *pNei = 0;
    //bool own_changed = false;
    //bool nei_changed = false;
    E_Int ao = -1, bo = -1, co = -1, an = -1, bn = -1, cn = -1;
    E_Int PHi_to_mod = -1;

    while (!stk.empty()) {
        PGi = stk.top();
        stk.pop();

        own = hmesh._F2E(0, PGi);
        nei = hmesh._F2E(1, PGi);

        if (nei == IDX_NONE) continue; // for now...
        if (!hmesh._PHtree.is_enabled(nei)) continue; // fow now i said...

        // make own left
        pOwn = hmesh._ng.PHs.get_facets_ptr(own);
        while (pOwn[3] != PGi+1) K_CONNECT::IdTool::right_shift<6>(pOwn, 1);
        K_MESH::Hexahedron::reorder_pgs_right(hmesh._ng, hmesh._F2E, own, 0);
        compute_canon_info_right(own, 0);
        Hexa_adap_compute(own, adap_incr);

        // make nei right
        pNei = hmesh._ng.PHs.get_facets_ptr(nei);
        while (pNei[2] != PGi+1) K_CONNECT::IdTool::right_shift<6>(pNei, 1);
        K_MESH::Hexahedron::reorder_pgs_left(hmesh._ng, hmesh._F2E, nei, 0);
        compute_canon_info_left(nei, 0);
        Hexa_adap_compute(nei, adap_incr);

        // store face_adap_incr for own and nei
        face_incr_t own_faces_incr[6];
        face_incr_t nei_faces_incr[6];

        for (E_Int i = 0; i < 6; i++) {
            own_faces_incr[i] = adap_incr.face_adap_incr[pOwn[i]-1];
            nei_faces_incr[i] = adap_incr.face_adap_incr[pNei[i]-1];
        }

        //own_changed = false;
        //nei_changed = false;

        //start_fix:
        ao = adap_incr.cell_adap_incr[own].n[0];
        bo = adap_incr.cell_adap_incr[own].n[1];
        co = adap_incr.cell_adap_incr[own].n[2];
        an = adap_incr.cell_adap_incr[nei].n[0];
        bn = adap_incr.cell_adap_incr[nei].n[1];
        cn = adap_incr.cell_adap_incr[nei].n[2];

        PHi_to_mod = -1;

        // fix alphas
        if (::abs(ao-an) <= 1) goto fix_beta;
        PHi_to_mod = ao > an ? nei : own;
        adap_incr.cell_adap_incr[PHi_to_mod].n[0] = std::max(ao, an) - 1;

        // fix betas
        fix_beta:
        if (::abs(bo-bn) <= 1) goto fix_gamma;
        PHi_to_mod = bo > bn ? nei : own;
        adap_incr.cell_adap_incr[PHi_to_mod].n[1] = std::max(bo, bn) - 1;

        // fix gammas
        fix_gamma:
        if (::abs(co-cn) <= 1) goto FIX;
        PHi_to_mod = co > cn ? nei : own;
        adap_incr.cell_adap_incr[PHi_to_mod].n[2] = std::max(co, cn) - 1;

        FIX:
        fix_beta_gamma(own, nei, adap_incr);

        // update own's and nei's faces
        adap_incr.face_adap_incr[PGi].n[0] = 0;
        adap_incr.face_adap_incr[PGi].n[1] = 0;
        fix_faces(own, adap_incr, hmesh);
        fix_faces(nei, adap_incr, hmesh);

        if (is_HX27(own, adap_incr)) {
            auto cell_incr = adap_incr.cell_adap_incr[own];
            E_Int count = 3;
            while (count == 3) {
                --cell_incr;
                count = 0;
                if (cell_incr.n[0]) count++;
                if (cell_incr.n[1]) count++;
                if (cell_incr.n[2]) count++;
            }

            if (count == 2) {
                // HX18
                
            }
        }

        // enforce consistancy between faces' increments and cells' increments
        // for each cell, check for existance of incomplete cutting patterns, and enforce them
        //if (!enforce_consistancy(own, hmesh, adap_incr) || !enforce_consistancy(nei, hmesh, adap_incr))
        //    goto start_fix;

        //check_consistancy(own, hmesh, adap_incr);
        //check_consistancy(nei, hmesh, adap_incr);

        // add neighbours to stack if shared face has changed
        // own
        for (E_Int i = 0; i < 6; i++) {
            if (i == 3) continue;
            E_Int face = pOwn[i]-1;
            if (adap_incr.face_adap_incr[face] != own_faces_incr[i])
                stk.push(face);
        }

        // nei
        for (E_Int i = 0; i < 6; i++) {
            if (i == 2) continue;
            E_Int face = pNei[i]-1;
            if (adap_incr.face_adap_incr[face] != nei_faces_incr[i])
                stk.push(face);
        }

        //std::cout << "stack size: " << stk.size() << std::endl;
    }

    for (PGi = 0; PGi < hmesh._ng.PGs.size(); PGi++) {
        if (hmesh._F2E(1,PGi) == IDX_NONE) {
            own = hmesh._F2E(0, PGi);
            pOwn = hmesh._ng.PHs.get_facets_ptr(own);
            while (pOwn[3] != PGi+1) K_CONNECT::IdTool::right_shift<6>(pOwn, 1);
            K_MESH::Hexahedron::reorder_pgs_right(hmesh._ng, hmesh._F2E, own, 0);
            compute_canon_info_right(own, 0);
            Hexa_adap_compute(own, adap_incr);

            auto& f_incr = adap_incr.face_adap_incr[PGi];
            f_incr.n[0] = 0; f_incr.n[1] = 0;
            E_Int beta = adap_incr.cell_adap_incr[own].n[1];
            E_Int gamma = adap_incr.cell_adap_incr[own].n[2];
            const auto& swap = _canon_info[own];
            if (beta) {
                if (swap[3]) f_incr.n[1] = beta;
                else         f_incr.n[0] = beta;
            }

            if (gamma) {
                if (swap[3]) f_incr.n[0] = gamma;
                else         f_incr.n[1] = gamma;
            }

            Hexa_adap_compute(own, adap_incr);
            assert(beta == adap_incr.cell_adap_incr[own].n[1]);
            assert(gamma == adap_incr.cell_adap_incr[own].n[2]);
        }
    }

    // reorder stuff
    E_Int alpha, beta, gamma;
    E_Int *pF;
    //E_Int cell_to_ref = 0;
    for (E_Int PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (adap_incr.cell_adap_incr[PHi] == 0) continue;

        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);

        assert(adap_incr.cell_adap_incr[PHi] != 0);

        alpha = adap_incr.cell_adap_incr[PHi].n[0];
        beta = adap_incr.cell_adap_incr[PHi].n[1];
        gamma = adap_incr.cell_adap_incr[PHi].n[2];

        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        //cell_to_ref++;
        
        if (is_HX12(PHi, adap_incr)) {
            if (alpha) {
                assert(!is_FXY(pF[0]-1, adap_incr));
                assert(!is_FXY(pF[1]-1, adap_incr));
                assert(!is_FXY(pF[4]-1, adap_incr));
                assert(!is_FXY(pF[5]-1, adap_incr));
            } else if (beta) {
                assert(!is_FXY(pF[0]-1, adap_incr));
                assert(!is_FXY(pF[1]-1, adap_incr));
                assert(!is_FXY(pF[2]-1, adap_incr));
                assert(!is_FXY(pF[3]-1, adap_incr));
            } else if (gamma) {
                K_CONNECT::IdTool::right_shift<6>(pF, 2);
                K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
                compute_canon_info_bottom(PHi, 0);
                Hexa_adap_compute(PHi, adap_incr);
                assert(is_HX12(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
                alpha = adap_incr.cell_adap_incr[PHi].n[0];
                beta = adap_incr.cell_adap_incr[PHi].n[1];
                if (alpha) {
                    assert(!is_FXY(pF[0]-1, adap_incr));
                    assert(!is_FXY(pF[1]-1, adap_incr));
                    assert(!is_FXY(pF[4]-1, adap_incr));
                    assert(!is_FXY(pF[5]-1, adap_incr));
                } else if (beta) {
                    assert(!is_FXY(pF[0]-1, adap_incr));
                    assert(!is_FXY(pF[1]-1, adap_incr));
                    assert(!is_FXY(pF[2]-1, adap_incr));
                    assert(!is_FXY(pF[3]-1, adap_incr));
                } else {
                    assert(false);
                }
            }
        } else if (is_HX18(PHi, adap_incr)) {
            if ((beta && gamma) || (alpha && gamma)) {
                while (!is_FXY(pF[0]-1, adap_incr)) K_CONNECT::IdTool::right_shift<6>(pF, 1);
                K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
                compute_canon_info_bottom(PHi, 0);
                Hexa_adap_compute(PHi, adap_incr);
                assert(is_HX18(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
            }
            assert(is_FXY(pF[0]-1, adap_incr));
            assert(is_FXY(pF[1]-1, adap_incr));
            assert(!is_FXY(pF[2]-1, adap_incr));
            assert(!is_FXY(pF[3]-1, adap_incr));
            assert(!is_FXY(pF[4]-1, adap_incr));
            assert(!is_FXY(pF[5]-1, adap_incr));
        }
    }
}

template <typename mesh_t>
void metric_sensor<mesh_t>::check_consistancy(E_Int PHi, mesh_t& hmesh, output_t& adap_incr)
{
    E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);
    const auto& swap = _canon_info[PHi];

    auto bot_incr(adap_incr.face_adap_incr[pF[0]-1]);
    auto top_incr(adap_incr.face_adap_incr[pF[1]-1]);
    auto lft_incr(adap_incr.face_adap_incr[pF[2]-1]);
    auto rgt_incr(adap_incr.face_adap_incr[pF[3]-1]);
    auto fro_incr(adap_incr.face_adap_incr[pF[4]-1]);
    auto bck_incr(adap_incr.face_adap_incr[pF[5]-1]);

    if (swap[0]) std::swap(bot_incr.n[0], bot_incr.n[1]);
    if (swap[1]) std::swap(top_incr.n[0], top_incr.n[1]);
    if (swap[2]) std::swap(lft_incr.n[0], lft_incr.n[1]);
    if (swap[3]) std::swap(rgt_incr.n[0], rgt_incr.n[1]);
    if (swap[4]) std::swap(fro_incr.n[0], fro_incr.n[1]);
    if (swap[5]) std::swap(bck_incr.n[0], bck_incr.n[1]);

    E_Int a = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int b = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int c = adap_incr.cell_adap_incr[PHi].n[2];

    assert(a == bot_incr.n[0]);
    assert(a == top_incr.n[0]);
    assert(a == fro_incr.n[0]);
    assert(a == bck_incr.n[0]);

    assert(b == bot_incr.n[1]);
    assert(b == top_incr.n[1]);
    assert(b == lft_incr.n[0]);
    assert(b == rgt_incr.n[0]);

    assert(c == lft_incr.n[1]);
    assert(c == rgt_incr.n[1]);
    assert(c == fro_incr.n[1]);
    assert(c == bck_incr.n[1]);
}

/*
template <typename mesh_t>
void metric_sensor<mesh_t>::metric_fix(mesh_t& hmesh, output_t& adap_incr)
{
    using cell_incr_t = typename output_t::cell_incr_t;
    using face_incr_t = typename output_t::face_incr_t;

    std::stack<E_Int> stk;
    E_Int ncells = hmesh._ng.PHs.size();

    for (E_Int PHi = 0; PHi < ncells; PHi++) {
        if (adap_incr.cell_adap_incr[PHi] != 0)
            stk.push(PHi);
    }

    while (!stk.empty()) {
        E_Int PHi = stk.top();
        stk.pop();

        // put PHi in canon config
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);

        E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        // prepare and store old_cell_incr
        E_Int BOT = prepare_bot_nei(PHi, hmesh, pF[0]-1, _start_nodes[PHi][0], adap_incr);
        E_Int TOP = prepare_top_nei(PHi, hmesh, pF[1]-1, _start_nodes[PHi][1], adap_incr);
        E_Int LFT = prepare_lft_nei(PHi, hmesh, pF[2]-1, _start_nodes[PHi][2], adap_incr);
        E_Int RGT = prepare_rgt_nei(PHi, hmesh, pF[3]-1, _start_nodes[PHi][3], adap_incr);
        E_Int FRO = prepare_fro_nei(PHi, hmesh, pF[4]-1, _start_nodes[PHi][4], adap_incr);
        E_Int BCK = prepare_bck_nei(PHi, hmesh, pF[5]-1, _start_nodes[PHi][5], adap_incr);

        cell_incr_t old_incr_BOT(cell_incr_t(0));
        cell_incr_t old_incr_TOP(cell_incr_t(0));
        cell_incr_t old_incr_LFT(cell_incr_t(0));
        cell_incr_t old_incr_RGT(cell_incr_t(0));
        cell_incr_t old_incr_FRO(cell_incr_t(0));
        cell_incr_t old_incr_BCK(cell_incr_t(0));

        if (BOT != IDX_NONE) old_incr_BOT = adap_incr.cell_adap_incr[BOT];
        if (TOP != IDX_NONE) old_incr_TOP = adap_incr.cell_adap_incr[TOP];
        if (LFT != IDX_NONE) old_incr_LFT = adap_incr.cell_adap_incr[LFT];
        if (RGT != IDX_NONE) old_incr_RGT = adap_incr.cell_adap_incr[RGT];
        if (FRO != IDX_NONE) old_incr_FRO = adap_incr.cell_adap_incr[FRO];
        if (BCK != IDX_NONE) old_incr_BCK = adap_incr.cell_adap_incr[BCK];

        bool changed = false;
        do {
            changed = false;
            changed |= FIX_BOT_TOP(PHi, BOT, adap_incr);
            changed |= FIX_BOT_TOP(PHi, TOP, adap_incr);
            changed |= FIX_LFT_RGT(PHi, LFT, adap_incr);
            changed |= FIX_LFT_RGT(PHi, RGT, adap_incr);
            changed |= FIX_FRO_BCK(PHi, FRO, adap_incr);
            changed |= FIX_FRO_BCK(PHi, BCK, adap_incr);
        } while (changed);

        // reset PHi's faces
        for (E_Int i = 0; i < 6; i++)
            adap_incr.face_adap_incr[pF[i]-1] = 0;
        fix_faces(PHi, adap_incr, hmesh);
        fix_faces(BOT, adap_incr, hmesh);
        fix_faces(TOP, adap_incr, hmesh);
        fix_faces(LFT, adap_incr, hmesh);
        fix_faces(RGT, adap_incr, hmesh);
        fix_faces(FRO, adap_incr, hmesh);
        fix_faces(BCK, adap_incr, hmesh);

        if ((BOT != IDX_NONE) && (old_incr_BOT != adap_incr.cell_adap_incr[BOT])) stk.push(BOT);
        if ((TOP != IDX_NONE) && (old_incr_TOP != adap_incr.cell_adap_incr[TOP])) stk.push(TOP);
        if ((LFT != IDX_NONE) && (old_incr_LFT != adap_incr.cell_adap_incr[LFT])) stk.push(LFT);
        if ((RGT != IDX_NONE) && (old_incr_RGT != adap_incr.cell_adap_incr[RGT])) stk.push(RGT);
        if ((FRO != IDX_NONE) && (old_incr_FRO != adap_incr.cell_adap_incr[FRO])) stk.push(FRO);
        if ((BCK != IDX_NONE) && (old_incr_BCK != adap_incr.cell_adap_incr[BCK])) stk.push(BCK);
    }

    E_Int alpha, beta, gamma;
    E_Int *pF;
    //E_Int cell_to_ref = 0;
    for (E_Int PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (adap_incr.cell_adap_incr[PHi] == 0) continue;

        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);

        assert(adap_incr.cell_adap_incr[PHi] != 0);

        alpha = adap_incr.cell_adap_incr[PHi].n[0];
        beta = adap_incr.cell_adap_incr[PHi].n[1];
        gamma = adap_incr.cell_adap_incr[PHi].n[2];

        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        //cell_to_ref++;
        
        if (is_HX12(PHi, adap_incr)) {
            if (alpha) {
                assert(!is_FXY(pF[0]-1, adap_incr));
                assert(!is_FXY(pF[1]-1, adap_incr));
                assert(!is_FXY(pF[4]-1, adap_incr));
                assert(!is_FXY(pF[5]-1, adap_incr));
            } else if (beta) {
                assert(!is_FXY(pF[0]-1, adap_incr));
                assert(!is_FXY(pF[1]-1, adap_incr));
                assert(!is_FXY(pF[2]-1, adap_incr));
                assert(!is_FXY(pF[3]-1, adap_incr));
            } else if (gamma) {
                K_CONNECT::IdTool::right_shift<6>(pF, 2);
                K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
                compute_canon_info_bottom(PHi, 0);
                Hexa_adap_compute(PHi, adap_incr);
                assert(is_HX12(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
                alpha = adap_incr.cell_adap_incr[PHi].n[0];
                beta = adap_incr.cell_adap_incr[PHi].n[1];
                if (alpha) {
                    assert(!is_FXY(pF[0]-1, adap_incr));
                    assert(!is_FXY(pF[1]-1, adap_incr));
                    assert(!is_FXY(pF[4]-1, adap_incr));
                    assert(!is_FXY(pF[5]-1, adap_incr));
                } else if (beta) {
                    assert(!is_FXY(pF[0]-1, adap_incr));
                    assert(!is_FXY(pF[1]-1, adap_incr));
                    assert(!is_FXY(pF[2]-1, adap_incr));
                    assert(!is_FXY(pF[3]-1, adap_incr));
                } else {
                    assert(false);
                }
            }
        } else if (is_HX18(PHi, adap_incr)) {
            if ((beta && gamma) || (alpha && gamma)) {
                while (!is_FXY(pF[0]-1, adap_incr)) K_CONNECT::IdTool::right_shift<6>(pF, 1);
                K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);
                compute_canon_info_bottom(PHi, 0);
                Hexa_adap_compute(PHi, adap_incr);
                assert(is_HX18(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
            }
            assert(is_FXY(pF[0]-1, adap_incr));
            assert(is_FXY(pF[1]-1, adap_incr));
            assert(!is_FXY(pF[2]-1, adap_incr));
            assert(!is_FXY(pF[3]-1, adap_incr));
            assert(!is_FXY(pF[4]-1, adap_incr));
            assert(!is_FXY(pF[5]-1, adap_incr));
        }
    }

    // loop on cells and check that face cuts match cell_adap_incr for HX12 and HX18
    for (E_Int PHi = 0; PHi < ncells; PHi++) {
        const auto& swap = _canon_info[PHi];
        const E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);
        E_Int alpha = adap_incr.cell_adap_incr[PHi].n[0];
        E_Int beta = adap_incr.cell_adap_incr[PHi].n[1];
        E_Int gamma = adap_incr.cell_adap_incr[PHi].n[2];
        E_Int BOT = pF[0]-1;
        E_Int TOP = pF[1]-1;
        E_Int LFT = pF[2]-1;
        E_Int RGT = pF[3]-1;
        E_Int FRO = pF[4]-1;
        E_Int BCK = pF[5]-1;
        const auto& bot_incr = adap_incr.face_adap_incr[BOT];
        const auto& top_incr = adap_incr.face_adap_incr[TOP];
        const auto& lft_incr = adap_incr.face_adap_incr[LFT];
        const auto& rgt_incr = adap_incr.face_adap_incr[RGT];
        const auto& fro_incr = adap_incr.face_adap_incr[FRO];
        const auto& bck_incr = adap_incr.face_adap_incr[BCK];
        if (is_HX12(PHi, adap_incr)) {
            assert(gamma == 0);
            if (alpha) {
                assert(!is_FXY(BOT, adap_incr));
                if (swap[0]) assert(bot_incr.n[1]);
                else assert(bot_incr.n[0]);

                assert(!is_FXY(TOP, adap_incr));
                if (swap[1]) assert(top_incr.n[1]);
                else assert(top_incr.n[0]);

                assert(!is_FXY(FRO, adap_incr));
                if (swap[4]) assert(fro_incr.n[1]);
                else assert(fro_incr.n[0]);

                assert(!is_FXY(BCK, adap_incr));
                if (swap[5]) assert(bck_incr.n[1]);
                else assert(bck_incr.n[0]);
            } else if (beta) {
                assert(!is_FXY(BOT, adap_incr));
                if (swap[0]) assert(bot_incr.n[0]);
                else assert(bot_incr.n[1]);

                assert(!is_FXY(TOP, adap_incr));
                if (swap[1]) assert(top_incr.n[0]);
                else assert(top_incr.n[1]);

                assert(!is_FXY(LFT, adap_incr));
                if (swap[2]) assert(lft_incr.n[1]);
                else assert(lft_incr.n[0]);

                assert(!is_FXY(RGT, adap_incr));
                if (swap[3]) assert(rgt_incr.n[1]);
                else assert(rgt_incr.n[0]);
            }
        } else if (is_HX18(PHi, adap_incr)) {
            assert(gamma == 0);
            assert(!is_FXY(LFT, adap_incr));
            if (swap[2]) assert(lft_incr.n[1]);
            else assert(lft_incr.n[0]);

            assert(!is_FXY(RGT, adap_incr));
            if (swap[3]) assert(rgt_incr.n[1]);
            else assert(rgt_incr.n[0]);

            assert(!is_FXY(FRO, adap_incr));
            if (swap[4]) assert(fro_incr.n[1]);
            else assert(fro_incr.n[0]);

            assert(!is_FXY(BCK, adap_incr));
            if (swap[5]) assert(bck_incr.n[1]);
            else assert(bck_incr.n[0]);
        }
    }
} // end function
*/

template <typename mesh_t>
void metric_sensor<mesh_t>::fix_faces(E_Int PHi, output_t& adap_incr, mesh_t& hmesh)
{
    if (PHi == IDX_NONE) return;

    const E_Int *faces = hmesh._ng.PHs.get_facets_ptr(PHi);
    
    E_Int bot = faces[0]-1;
    E_Int top = faces[1]-1;
    E_Int lft = faces[2]-1;
    E_Int rgt = faces[3]-1;
    E_Int fro = faces[4]-1;
    E_Int bck = faces[5]-1;

    auto& bot_incr = adap_incr.face_adap_incr[bot];
    auto& top_incr = adap_incr.face_adap_incr[top];
    auto& lft_incr = adap_incr.face_adap_incr[lft];
    auto& rgt_incr = adap_incr.face_adap_incr[rgt];
    auto& fro_incr = adap_incr.face_adap_incr[fro];
    auto& bck_incr = adap_incr.face_adap_incr[bck];

    const auto& swap = _canon_info[PHi];

    E_Int alpha = adap_incr.cell_adap_incr[PHi].n[0];
    E_Int beta = adap_incr.cell_adap_incr[PHi].n[1];
    E_Int gamma = adap_incr.cell_adap_incr[PHi].n[2];

    // alpha -> botX, topX, froX, bckX
    if (swap[0]) bot_incr.n[1] = std::max(bot_incr.n[1], alpha);
    else         bot_incr.n[0] = std::max(bot_incr.n[0], alpha);
    if (swap[1]) top_incr.n[1] = std::max(top_incr.n[1], alpha);
    else         top_incr.n[0] = std::max(top_incr.n[0], alpha);
    if (swap[4]) fro_incr.n[1] = std::max(fro_incr.n[1], alpha);
    else         fro_incr.n[0] = std::max(fro_incr.n[0], alpha);
    if (swap[5]) bck_incr.n[1] = std::max(bck_incr.n[1], alpha);
    else         bck_incr.n[0] = std::max(bck_incr.n[0], alpha);

    // beta -> botY, topY, lftX, rgtX
    if (swap[0]) bot_incr.n[0] = std::max(bot_incr.n[0], beta);
    else         bot_incr.n[1] = std::max(bot_incr.n[1], beta);
    if (swap[1]) top_incr.n[0] = std::max(top_incr.n[0], beta);
    else         top_incr.n[1] = std::max(top_incr.n[1], beta);
    if (swap[2]) lft_incr.n[1] = std::max(lft_incr.n[1], beta);
    else         lft_incr.n[0] = std::max(lft_incr.n[0], beta);
    if (swap[3]) rgt_incr.n[1] = std::max(rgt_incr.n[1], beta);
    else         rgt_incr.n[0] = std::max(rgt_incr.n[0], beta);

    // gamma -> lftY, rgtY, froY, bckY
    if (swap[2]) lft_incr.n[0] = std::max(lft_incr.n[0], gamma);
    else         lft_incr.n[1] = std::max(lft_incr.n[1], gamma);
    if (swap[3]) rgt_incr.n[0] = std::max(rgt_incr.n[0], gamma);
    else         rgt_incr.n[1] = std::max(rgt_incr.n[1], gamma);
    if (swap[4]) fro_incr.n[0] = std::max(fro_incr.n[0], gamma);
    else         fro_incr.n[1] = std::max(fro_incr.n[1], gamma);
    if (swap[5]) bck_incr.n[0] = std::max(bck_incr.n[0], gamma);
    else         bck_incr.n[1] = std::max(bck_incr.n[1], gamma);
}

template <typename mesh_t>
void metric_sensor<mesh_t>::metric_fix_2(mesh_t& hmesh, output_t& adap_incr)
{
    using cell_incr_t = typename output_t::cell_incr_t;
    //using face_incr_t = typename output_t::face_incr_t;

    cell_incr_t c_incr_1, c_incr_0;
    E_Int count;
    E_Int *pF = 0;
    E_Int alpha, beta, gamma;
    std::stack<E_Int> stk;
    E_Int PGi, PHi;

    for (PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (adap_incr.cell_adap_incr[PHi] == 0) continue;
        stk.push(PHi);
    }

    while (!stk.empty()) {
        PHi = stk.top();
        stk.pop();

        count = 0;
        for (E_Int i = 0; i < 3; i++) {
            if (adap_incr.cell_adap_incr[PHi].n[i])
                count++;
        }

        if (count == 0 || count == 3) continue;

        alpha = adap_incr.cell_adap_incr[PHi].n[0];
        beta = adap_incr.cell_adap_incr[PHi].n[1];
        gamma = adap_incr.cell_adap_incr[PHi].n[2];

        const auto& swap = _canon_info[PHi];
        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        for (int i = 0; i < 6; i++) {
            PGi = pF[i]-1;
            if (!is_FXY(PGi, adap_incr)) continue;

            E_Int nei = NEIGHBOR(PHi, hmesh._F2E, PGi);
            bool ok_nei = (nei != IDX_NONE) && hmesh._PHtree.is_enabled(PHi);

            auto& PG_incr = adap_incr.face_adap_incr[PGi];

            if (count == 1) {
                if (alpha) {
                    if (i != 2 && i != 3) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (beta) {
                    if (i == 0 || i == 1) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    } else if (i == 2 || i == 3) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (gamma) {
                    if (i != 0 && i != 1) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                }
            } else if (count == 2) {
                if (alpha && beta) {
                    if (i != 0 && i != 1) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (beta && gamma) {
                    if (i != 2 && i != 3) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                } else if (gamma && alpha) {
                    if (i == 0 || i == 1) {
                        if (swap[i]) PG_incr.n[0] = 0;
                        else PG_incr.n[1] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    } else if (i == 2 || i == 3) {
                        if (swap[i]) PG_incr.n[1] = 0;
                        else PG_incr.n[0] = 0;
                        if (ok_nei) {
                            Hexa_adap_compute(nei, adap_incr);
                            if (adap_incr.cell_adap_incr[nei] != 0) stk.push(nei);
                        }
                    }
                }
            }
        }
    }

    

    // reorder HX12 in gamma to alpha
    for (PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (!is_HX12(PHi, adap_incr)) continue;
        gamma = adap_incr.cell_adap_incr[PHi].n[2];
        if (gamma == 0) continue;

        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        K_CONNECT::IdTool::right_shift<6>(pF, 2);
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);

        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);
        assert(is_HX12(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
    }

    // reorder HX18
    for (PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (!hmesh._PHtree.is_enabled(PHi)) continue;
        if (!is_HX18(PHi, adap_incr)) continue;

        pF = hmesh._ng.PHs.get_facets_ptr(PHi);

        E_Int i0 = 0;
        while (!is_FXY(pF[0]-1, adap_incr)) {
            K_CONNECT::IdTool::right_shift<6>(pF, 1);
            i0++;
            if (i0 == 6) {
                std::cout << "oups\n";
                break;
            }
        }
            
        K_MESH::Hexahedron::reorder_pgs(hmesh._ng, hmesh._F2E, PHi);

        compute_canon_info_bottom(PHi, 0);
        Hexa_adap_compute(PHi, adap_incr);
        assert(is_HX18(PHi, adap_incr) && adap_incr.cell_adap_incr[PHi].n[2] == 0);
    }

    for (E_Int PHi = 0; PHi < hmesh._ng.PHs.size(); PHi++) {
        if (adap_incr.cell_adap_incr[PHi] == 0) continue;
        const E_Int *pF = hmesh._ng.PHs.get_facets_ptr(PHi);
        E_Int alpha = adap_incr.cell_adap_incr[PHi].n[0];
        E_Int beta = adap_incr.cell_adap_incr[PHi].n[1];
        E_Int gamma = adap_incr.cell_adap_incr[PHi].n[2];

        E_Int BOT = pF[0]-1;
        E_Int TOP = pF[1]-1;
        E_Int LFT = pF[2]-1;
        E_Int RGT = pF[3]-1;
        E_Int FRO = pF[4]-1;
        E_Int BCK = pF[5]-1;

        if (is_HX12(PHi, adap_incr)) {
            if (alpha) {
                assert(!is_FXY(BOT, adap_incr));
                assert(!is_FXY(TOP, adap_incr));
                assert(!is_FXY(FRO, adap_incr));
                assert(!is_FXY(BCK, adap_incr));
            } else if (beta) {
                assert(!is_FXY(BOT, adap_incr));
                assert(!is_FXY(TOP, adap_incr));
                assert(!is_FXY(LFT, adap_incr));
                assert(!is_FXY(RGT, adap_incr));
            } else if (gamma) {
                assert(!is_FXY(LFT, adap_incr));
                assert(!is_FXY(RGT, adap_incr));
                assert(!is_FXY(FRO, adap_incr));
                assert(!is_FXY(BCK, adap_incr));
            }
        } else if (is_HX18(PHi, adap_incr)) {
             if (alpha && beta) {
                assert(is_FXY(BOT, adap_incr));
                assert(is_FXY(TOP, adap_incr));
                assert(!is_FXY(FRO, adap_incr));
                assert(!is_FXY(BCK, adap_incr));
                assert(!is_FXY(LFT, adap_incr));
                assert(!is_FXY(RGT, adap_incr));
            } else if (alpha && gamma) {
                assert(is_FXY(FRO, adap_incr));
                assert(is_FXY(BCK, adap_incr));
                assert(!is_FXY(LFT, adap_incr));
                assert(!is_FXY(RGT, adap_incr));
                assert(!is_FXY(TOP, adap_incr));
                assert(!is_FXY(BOT, adap_incr));
            } else if (beta && gamma) {
                assert(is_FXY(LFT, adap_incr));
                assert(is_FXY(RGT, adap_incr));
                assert(!is_FXY(FRO, adap_incr));
                assert(!is_FXY(BCK, adap_incr));
                assert(!is_FXY(BOT, adap_incr));
                assert(!is_FXY(TOP, adap_incr));
            }
        }
    }

    std::cout << "adap_incr fixed\n";
}

/// 
template <typename mesh_t>
E_Int metric_sensor<mesh_t>::assign_data(const sensor_input_t& data)
{
  E_Int npts = parent_t::_hmesh._crd.cols();
  if (npts != data.size()) exit(1);

  parent_t::_data = data;

  _single_pass_done = false;//reinit to enable

  return 0;
}

template <typename mesh_t>
void metric_sensor<mesh_t>::Q4_adap_compute(E_Int PGi, output_t& adap_incr)
{
    //int h[4];
    E_Float h[4];
    E_Int n1, n2;
    auto& metric = parent_t::_data;
    const E_Int *pN = parent_t::_hmesh._ng.PGs.get_facets_ptr(PGi);

    for (int i = 0; i < 4; i++) {
        n1 = pN[i]; n2 = pN[(i+1)%4];
        h[i] = metric.lengthEval(n1-1, metric[n1-1], n2-1, metric[n2-1]);
    }

    auto f1 = log2(std::min(h[0], h[2]));
    auto f2 = log2(std::min(h[1], h[3]));

    adap_incr.face_adap_incr[PGi].n[0] = round(f1);
    adap_incr.face_adap_incr[PGi].n[1] = round(f2);
}

template <typename mesh_t>
void metric_sensor<mesh_t>::Hexa_adap_compute(E_Int PHi, output_t& adap_incr)
{
    if (PHi == IDX_NONE) return;
    if (!parent_t::_hmesh._PHtree.is_enabled(PHi)) return;

    using face_incr_t = typename output_t::face_incr_t;

    const auto& ng = parent_t::_hmesh._ng;

    face_incr_t f_incr[6];
    E_Int PGi;
    const E_Int *pF = ng.PHs.get_facets_ptr(PHi);
    const auto& swap = _canon_info[PHi];

    for (int i = 0; i < 6; i++) {
        PGi = pF[i] - 1;
        if (swap[i] > 0) {
            f_incr[i].n[0] = adap_incr.face_adap_incr[PGi].n[1];
            f_incr[i].n[1] = adap_incr.face_adap_incr[PGi].n[0];
        } else {
            f_incr[i].n[0] = adap_incr.face_adap_incr[PGi].n[0];
            f_incr[i].n[1] = adap_incr.face_adap_incr[PGi].n[1];
        }
    }

    using cell_incr_t = typename output_t::cell_incr_t;

    cell_incr_t cell_incr(0);

    E_Int incrX[4] = {f_incr[0].n[0], f_incr[4].n[0], f_incr[1].n[0], f_incr[5].n[0]};
    cell_incr.n[0] = std::max(E_Int(0),*std::min_element(incrX, incrX + 4));

    E_Int incrY[4] = {f_incr[0].n[1], f_incr[2].n[0], f_incr[1].n[1], f_incr[3].n[0]};
    cell_incr.n[1] = std::max(E_Int(0),*std::min_element(incrY, incrY + 4));

    E_Int incrZ[4] = {f_incr[2].n[1], f_incr[4].n[1], f_incr[3].n[1], f_incr[5].n[1]};
    cell_incr.n[2] = std::max(E_Int(0),*std::min_element(incrZ, incrZ + 4));

    adap_incr.cell_adap_incr[PHi] = cell_incr;
}

///
template <typename mesh_t>
bool metric_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  using cell_incr_t = typename output_t::cell_incr_t; //int_tuple<3>
  using face_incr_t = typename output_t::face_incr_t; // int_tuple<2>
  
  // init face_adap_incr
  adap_incr.face_adap_incr.clear();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  // init cell_adap_incr
  adap_incr.cell_adap_incr.clear();
  E_Int nb_cells = parent_t::_hmesh._ng.PHs.size();
  adap_incr.cell_adap_incr.resize(nb_cells, cell_incr_t(0));

  _single_pass_done = true; // to exit after first outer loop (based on sensor update)

  E_Int PHi, PGi, stride;
  const auto& F2E = parent_t::_hmesh._F2E;
  E_Int nifaces = 0;
  for (E_Int PGi = 0; PGi < nb_faces; PGi++) {
      if (F2E(0, PGi) != IDX_NONE || F2E(1, PGi) != IDX_NONE) nifaces++;
  }
  std::cout << "nifaces: " << nifaces << std::endl;

  for (PHi = 0; PHi < nb_cells; PHi++) {
    if (!parent_t::_hmesh._PHtree.is_enabled(PHi)) continue;
    stride = parent_t::_hmesh._ng.PHs.stride(PHi);
    const E_Int *pF = parent_t::_hmesh._ng.PHs.get_facets_ptr(PHi);
    for (int j = 0; j < stride; j++) {
      PGi = pF[j] - 1;
      // only hexa for now
      Q4_adap_compute(PGi, adap_incr);
    }
    compute_canon_info_bottom(PHi, 0);
    Hexa_adap_compute(PHi, adap_incr);
  }
  std::cout << "adap_incr filled\n";

  metric_fix(parent_t::_hmesh, adap_incr);

  E_Int cmax = 0;
  for (size_t k = 0; (k < adap_incr.cell_adap_incr.size()) /*&& (cmax==0)*/; ++k)
    cmax = std::max(cmax, adap_incr.cmax(k));

  E_Int cmin = IDX_NONE; // max int32
  for (size_t k = 0; (k < adap_incr.cell_adap_incr.size()) && (cmin != 0); ++k)
    cmin = std::min(cmin, adap_incr.cmin(k));

  std::cout << "cell adapt incr min/max : " << cmin << "/" << cmax << std::endl;

  return (cmin != 0 || cmax != 0);
}

}
#endif
