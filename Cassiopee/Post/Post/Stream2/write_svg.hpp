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
#ifndef _SVG_HPP_
#define _SVG_HPP_
#include <cstdio>
#include <string>
#include <array>
#include "kcore.h"

namespace SVG
{
  struct point
  {
    E_Float _x, _y;
    point() = default;
    point(E_Float x, E_Float y) : _x(x), _y(y) {}
  };

  /** Ouvre un fichier en écriture pour écrire une image au format SVG */
  FILE* open_svg(const char* fileName);
  /** Rajoute une description pour le format SVG (un commentaire) */
  void write_description_svg(FILE* fich, const std::string& desc);
  /** Rajoute un point dans le fichier SVG (sous forme de petit rectangle) */
  void write_point_svg(FILE* fich, const point& p, const std::string&  col = "black" );
  /** Rajoute une ligne allant du point p1 au point p2 */
  void write_line_svg(FILE* fich, const point& p1, const point& p2, const std::string& color = "red");
  /** Rajoute un cercle creux de centre p et de rayon r */
  void write_circle_svg(FILE* fich, const point& p, E_Float r, int thick = 1, const std::string& color = "blue");
  /** Rajoute un cercle plein de centre p et de rayon r */
  void write_circle_svg(FILE* fich, const point& p, E_Float r, const std::string& fill_color, E_Float fill_opactity = 1.);
  /** Rajoute un rectangle creux ayant pour coins opposés p_min et p_max (coords min et max) */
  void write_rectangle_svg(FILE* fich, const point& p_min, const point& p_max, int thick = 1, const std::string& color = "red" );
  /** Rajoute un rectangle plein ayant pour coins opposés p_min et p_max (coords min et max) */
  void write_rectangle_svg(FILE* fich, const point& p_min, const point& p_max, const std::string& fill_color, E_Float fill_opacity = 1.);
  /** Rajoute un texte sur l'image avec un éventuel remplissage */
  void write_text_svg(FILE* fich, const point& origin, const std::string& texte, const std::string& fill_color = "black", E_Float fill_opacity = 1. );
  /** Affiche un triangle creux */
  void write_triangle_svg(FILE* fich, const point& p1, const point& p2, const point& p3, int thick = 1, const std::string& color = "red");
  /** Affiche un triangle plein */
  void write_triangle_svg(FILE* fich, const point& p1, const point& p2, const point& p3, const std::string& fill_color, E_Float fill_opacity = 1.);
  /** Referme un fichier svg avec les bonnes balises */
  void close_svg(FILE* fich);
}

#endif
