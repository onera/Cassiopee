#ifndef _CPLOT_SHADERS_ID_H_
#define _CPLOT_SHADERS_ID_H_

namespace shader
{
enum shaders_id
{
	None = 0,
	unidirectionnel_phong,
	bidirectionnel_phong,
	glass,
	chrome,
	metal,
	wood,
	marble,
	smoke,
	xray,
	iso_banded_colormap,
	granite,
	sphere_billboarding,
	monochrome_anaglyph,
	color_anaglyph,
	iso_continuous_colormap,
	brick,
	cloud,
	iso_granite,
	shadow_mapping,
	DOF,
	gooch,
	flat,
	billboard,
	iso_flat,
	iso_chrome,
	iso_glass,
	rgb_vector,
	iso_brick,
	iso_colored_lines,
	iso_xray,
	iso_gooch,
	iso_metal,
	velocity_line,
	velocity_tetra,
	velocity_triangle,
	velocity_uniform_streamline,
	textured_material,

	end_of_shaders_id
};
}
#endif