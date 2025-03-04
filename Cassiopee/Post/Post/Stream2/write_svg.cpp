/*
    Copyright 2013-2025 Onera.

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
#include "write_svg.hpp"

namespace
{
    inline double
    X(double x)
    {
        return (50. + x) / 10;
    }

    inline double
    Y(double y)
    {
        return (50. - y) / 10;
    }
}  // namespace

namespace SVG
{
    FILE*
    open_svg(const char* fileName)
    {
        FILE* fich = fopen(fileName, "w");
        fprintf(fich, "<?xml version=\"1.0\" standalone=\"no\"?>\n"),
            fprintf(fich,
                    "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\""
                    "\n \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        fprintf(fich, "<svg width=\"10cm\" height=\"10cm\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n");
        return fich;
    }

    void
    write_description_svg(FILE* fich, const std::string& desc)
    {
        fprintf(fich, "<desc>%s\n", desc.c_str());
        fprintf(fich, "</desc>\n");
    }

    void
    write_point_svg(FILE* fich, const point& p, const std::string& col)
    {
        fprintf(fich, "<rect x=\"%gcm\" y=\"%gcm\" width=\"0.1cm\" height=\"0.1cm\" fill=\"%s\"/>\n",
                (49.5 + p._x) / 10., (49.5 - p._y) / 10., col.c_str());
    }

    void
    write_line_svg(FILE* fich, const point& p1, const point& p2, const std::string& col)
    {
        fprintf(fich, "<line x1=\"%gcm\" y1=\"%gcm\" x2=\"%gcm\" y2=\"%gcm\" stroke=\"%s\"/>\n", X(p1._x), Y(p1._y),
                X(p2._x), Y(p2._y), col.c_str());
    }

    void
    write_circle_svg(FILE* fich, const point& p, double r, int thick, const std::string& color_border)
    {
        fprintf(fich, "<circle cx=\"%gcm\" cy=\"%gcm\" r=\"%gcm\" stroke=\"%s\" stroke-width=\"%d\" fill=\"none\" />",
                X(p._x), Y(p._y), r / 10., color_border.c_str(), thick);
    }

    void
    write_circle_svg(FILE* fich, const point& p, double r, const std::string& fill_color, double fill_opacity)
    {
        fprintf(fich,
                "<circle cx=\"%gcm\" cy=\"%gcm\" r=\"%gcm\" stroke-opacity=\"0\" fill=\"%s\" fill-opacity=\"%g\" />",
                X(p._x), Y(p._y), r / 10., fill_color.c_str(), fill_opacity);
    }

    void
    write_rectangle_svg(FILE* fich, const point& p_min, const point& p_max, int thick, const std::string& color)
    {
        fprintf(fich,
                "<rect x=\"%gcm\" y=\"%gcm\" width=\"%gcm\" height=\"%gcm\" stroke=\"%s\" stroke-width=\"%d\" "
                "fill=\"none\" />\n",
                X(p_min._x), Y(p_max._y), (p_max._x - p_min._x) / 10., (p_max._y - p_min._y) / 10., color.c_str(),
                thick);
    }

    void
    write_rectangle_svg(FILE* fich, const point& p_min, const point& p_max, const std::string& fill_color,
                        double fill_opacity)
    {
        fprintf(fich,
                "<rect x=\"%gcm\" y=\"%gcm\" width=\"%gcm\" height=\"%gcm\" stroke_opacity=\"0\" fill=\"%s\" "
                "fill-opacity=\"%g\" />\n",
                X(p_min._x), Y(p_max._y), (p_max._x - p_min._x) / 10., (p_max._y - p_min._y) / 10., fill_color.c_str(),
                fill_opacity);
    }

    void
    write_text_svg(FILE* fich, const point& p, const std::string& texte, const std::string& fill_color,
                   double fill_opacity)
    {
        fprintf(fich, "<text x=\"%gcm\" y=\"%gcm\" fill=\"%s\" fill-opacity=\"%g\">%s</text>\n", X(p._x), Y(p._y),
                fill_color.c_str(), fill_opacity, texte.c_str());
    }

    void
    write_triangle_svg(FILE* fich, const point& p1, const point& p2, const point& p3, int thick,
                       const std::string& color)
    {
        fprintf(fich, "<polygon points=\"%g,%g %g,%g %g,%g\" stroke=\"%s\" stroke-width=\"%d\" fill=\"none\" />\n",
                37.795 * X(p1._x), 37.795 * Y(p1._y), 37.795 * X(p2._x), 37.795 * Y(p2._y), 37.795 * X(p3._x),
                37.795 * Y(p3._y), color.c_str(), thick);
    }

    void
    write_triangle_svg(FILE* fich, const point& p1, const point& p2, const point& p3, const std::string& fill_color,
                       double fill_opacity)
    {
        fprintf(fich,
                "<g transform=\"scale(37.795)\">\n\t<polygon points=\"%g,%g %g,%g %g,%g\" fill=\"%s\" "
                "fill-opacity=\"%g\" />\n</g>\n",
                X(p1._x), Y(p1._y), X(p2._x), Y(p2._y), X(p3._x), Y(p3._y), fill_color.c_str(), fill_opacity);
    }


    void
    close_svg(FILE* fich)
    {
        fprintf(fich, "<!-- Show outline of canvas using 'rect' element -->\n");
        fprintf(fich,
                "<rect x=\".01cm\" y=\".01cm\" width=\"9.98cm\" height=\"9.98cm\""
                "\n fill=\"none\" stroke=\"blue\" stroke-width=\".02cm\" />\n");
        fprintf(fich, "</svg>\n");
        fclose(fich);
    }
}  // namespace SVG
