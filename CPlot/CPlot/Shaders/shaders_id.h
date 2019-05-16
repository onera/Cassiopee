#ifndef _CPLOT_SHADERS_ID_H_
#define _CPLOT_SHADERS_ID_H_

namespace shader
{
enum shaders_id
{
	None = 0,
	unidirectionnel_phong, // 1
	bidirectionnel_phong,  // 2
	glass, // 3
	chrome, // 4
	metal, // 5
	wood, // 6
	marble, // 7
	smoke, // 8
	xray, // 9
	iso_banded_colormap, // 10
	granite, // 11
	sphere_billboarding, // 12
	monochrome_anaglyph, // 13
	color_anaglyph, // 14
	iso_continuous_colormap, // 15
	brick, // 16
	cloud, // 17
	iso_granite, // 18
	shadow_mapping, // 19
	DOF, // 20
	gooch, // 21
	flat, // 22
	billboard, // 23
	iso_flat, // 24
	iso_chrome, // 25
	iso_glass, // 26
	vector_rgb, // 27
	iso_brick, // 28
	iso_colored_lines, // 29
	iso_xray, // 30
	iso_gooch, // 31
	iso_metal, // 32
	vector_line, // 33
	vector_arrow, // 34
	textured_material, // 35

	end_of_shaders_id // 36 !
};
}
#endif