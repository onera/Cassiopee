# Benchmark for generating parametric geometries

## parametric two circles

Read @tools/agents/AGENTS_user.md for generic instructions.
Write a python script that creates a 2D parametric geometry of two circles.
Choose the parameters to have about 3 parameters and that are key parameters to the shape.
Write all script and artifacs to directory @/tools/sandbox/circles.

semia-small: FAIL
semia-large: SUCCESS

## parametric rocket

Read @tools/agents/AGENTS_user.md for generic instructions.
Write a python script that creates a 2D parametric sketch of a rocket.
Choose the parameters to have about 5 parameters and that are key parameters to the shape.
Write all scripts and artifacs to directory @/tools/sandbox/rocket.

semia-small: 
semia-large: SUCESS

## Roms Letter C

Read @tools/agents/AGENTS_user.md for generic instructions.
Write a python script that creates a 2D parametric geometry 
of the letter "C".
Choose the parameters to have 3 parameters and that are key parameters to the shape.
The shape must be made of a single connected and closed sketch.
Write all scripts and artifacs to directory @/tools/sandbox/letterC.

semia-small: overflow context, infinite loop, FAIL
semia-large: SUCCESS (5 minutes)
qwen3.5: FAIL

## parametric pacman shape

Read @tools/agents/AGENTS_user.md for generic instructions.
Write a python script that create 
a 2D parametric sketch of a pacman shape.
Choose the parameters to have 4 parameters that are key parameters to the shape.
Write all scripts and artifacs to directory @/tools/sandbox/pacman.

semia-small: SUCCESS (4 minutes)
semia-large:

## parametric smiley shape

Read @tools/agents/AGENTS_user.md for generic instructions.
Write a python script that create 
a 2D parametric sketch of a
"smiley". The shape must be done of two circles
for the eyes and one spline for the mouth and one big circle around all.
Choose the parameters to have 4 parameters that are key parameters to the shape.
Write all script and artifacs to @/tools/sandbox/smiley.

semia-small: 
semia-large: SUCCESS (4 minutes)

## parametric cubes

Read @tools/agents/AGENTS_user.md for generic instructions.
Write a python script that create 
a 3D parametric surface made of the union of two cubes.
Choose the parameters to have 4 parameters that are key parameters to the shape.
Write all script and artifacs to directory @/tools/sandbox/cubes.

semia-small: 
semia-large: SUCCESS (6 minutes)
