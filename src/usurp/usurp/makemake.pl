#!/usr/bin/perl
#
# Usage: makemake {<program name> {<F90 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# This script recognizes the operating system and assumes that it comes
# with the vendor compilers.  Otherwise, must specify the 3rd party
# compilers.
# Modified by Birol Aktas <baktas@scientech.com> December 7, 1998
# SCIENTECH, Inc., Rockville, Maryland
#
#
# Modified 7/1/98 by Simon Smith to change defaults to debugging version
# and add *.M, *.mod
#
# Modified 5/2005 by David Boger for customization for USURP
#
$platform = `uname`;
chomp($platform);
if ($platform eq "CYGWIN_NT-4.0") {
  $o = "obj";
  $out = "/out:";
  $mod = "mod"
} else {
  $o = "o";
  $out = "-o ";
  $mod = "mod";
}
#
open(MAKEFILE, "> Makefile");
#
print MAKEFILE "PROG =\t$ARGV[0]\n\n";
#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <*.f90 *.F90 *.f *.F *.c>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.$o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Module listing
#
print MAKEFILE "MODS =\t";
foreach $file (<*.f90 *.F90>) {
  open(FILE, $file) || warn "Cannot open $file: $!\n";
  while (<FILE>) {
    if (/^\s*module\s+([^\s!]+)/i) {
      $filename = &toLower($1);
      if ($filename ne "procedure") {push(@mods, "$filename\.$mod")};
    }
  }
  close(FILE);
}
&PrintWords(8, 0, @mods);
print MAKEFILE "\n\n";

#
# Define common macros
#
print MAKEFILE "CMD =\tall\n\n";

if ($platform eq "OSF1") {
  $LIBS = "../Lib/libtracc.a ../Lib/libfpvm3.a ../Lib/libgpvm3.a ../Lib/libpvm3.a ../Lib/libdxml.a";
  $CC = "cc"; 
  $CFLAGS = "-i4 -r8 -ladebug -DDECOSF";
  $FC = "f77";
  $FFLAGS = "-g";
  $F90 = "f90";
  $F90FLAGS = "-c -g -cpp -fixed -i4 -r8 -ladebug -check_bounds -fpe2";
  $LDFLAGS = "";
} elsif ($platform eq "Linux") {
  $LIBS = "\$(TECIO_HOME) \$(FXDR_HOME)";
  $CC = "gcc";
  $CFLAGS = "-c -O2";
  $FC = "";
  $FFLAGS = "";
  $F90 = "pgf90";
  $FPPFLAGS = "";
  $F90FLAGS = "-c -O2 -Mstandard -Mbyteswapio -DUSE_TECPLOT -DMCLOCK";
  $LDFLAGS = "-O2";
} elsif ($platform eq "CYGWIN_NT-4.0") {
  $LIBS = "kernel32.lib ws2_32.lib user32.lib wsock32.lib cFiles.lib SuperLU.lib PvmLibs.lib ADVAPI32.LIB ";
  $CC = "cl";
  $CFLAGS = "/nologo -c";
  $VPATH = "Debug";
  $FC = "DF";
  $FFLAGS = "-g";
  $LINK = "link";
  $F90 = "f90";
  $F90FLAGS = "/nologo /compile_only /keep /warn:argument_checking /warn:declarations /nofree /check:bounds /Debug:full /module:'Debug' /object:'Debug' /browser:'Debug' /traceback /warn:nofileopt /warn:nouncalled /warn:nouninitialized /pdbfile:'Debug'";
  $LDFLAGS = "/nodefaultlib:'libcd' /libpath:'Debug' /libpath:'d:/work/advcode/cFiles/Debug' /libpath:'d:/work/advcode/SuperLU/Debug' /libpath:'d:/work/advcode/PvmLibs/Debug'";
} elsif ($platform eq "IRIX64") {
  $LIBS = "../Lib/libtracc.a ../Lib/libfpvm3.a ../Lib/libgpvm3.a ../Lib/libpvm3.a";
  $CC = "cc";
  $CFLAGS = "-g -DSGI -i4 -64 ";
  $FC = "f77";
  $FFLAGS = "-g";
  $F90 = "f90";
  $F90FLAGS = "-i4 -64 -c -cpp -g -woff878,1563,1582,399,1438 -DEBUG:div_check=3:trap_uninitialized=ON:verbose_runtime -fixedform -G0 -O0 -fullwarn";
  $LDFLAGS = "-i4 -64 -g -DEBUG:div_check=3:subscript_check=ON:trap_uninitialized=ON:verbose_runtime=ON -fixedform -check_bounds -G0 -O0";
} elsif ($platform eq "SunOS") {
  $LIBS = "";
  $CC = "cc";
  $CFLAGS = "-c -O";
  $FC = "f77";
  $FFLAGS = "-c -O";
  $F90 = "f90";
  $F90FLAGS = "-c -O -DCPU_TIME";
  $LDFLAGS = "";
} elsif ($platform eq "HPUX") {
  $LIBS = "";
  $CC = "cc"; 
  $CFLAGS = "-g -DHP -Aa";
  $FC = "f77";
  $FFLAGS = "-g";
  $F90 = "f90";
  $F90FLAGS = "-c -g";
  $LDFLAGS = "";
} elsif ($platform eq "AIX") {
  $LIBS = "/usr/cta/tecplot/10/lib/tecio64.a";
  $CC = "xlc"; 
  $CFLAGS = "-c -q64 -O2";
  $FC = "f77";
  $FFLAGS = "-c -q64 -O2";
  $F90 = "xlf90";
  $F90FLAGS = "-c -q64 -O2 -qlanglvl=90std -qsuffix=f=f90:cpp=F90 -qextname -WF,-DUSE_TECPLOT,-DMCLOCK";
  $LDFLAGS = "-q64 -O2 -bmaxdata:0x80000000";
# -qqfixed=72 -qmaxmem=16352 -qspillsize=16352"; for FFLAGS and F90FLAGS
# -qspillsize: the stack size of the argument lists
# -qmaxmem: controls the internal compiler memory 
} else {
  print "Platform unknown \n";
  exit
}
print MAKEFILE "LIBS = $LIBS\t\n\n";
print MAKEFILE "CC = $CC\n";
print MAKEFILE "CFLAGS = $CFLAGS\n";
print MAKEFILE "FC = $FC\n";
print MAKEFILE "FFLAGS = $FFLAGS\n";
print MAKEFILE "F90 = $F90\n";
print MAKEFILE "F90FLAGS = $F90FLAGS\n";
if ($platform eq "CYGWIN_NT-4.0") {
  print MAKEFILE "VPATH = $VPATH\n";
  print MAKEFILE "LINK = $LINK\n";
}
print MAKEFILE "LDFLAGS = $LDFLAGS\n\n";
#
# make
#
print MAKEFILE "include Make.sys\n\n";
print MAKEFILE "all: \${PROG}\n\n";
print MAKEFILE "\${PROG}: \${OBJS}\n";
#print MAKEFILE "\trm -f Debug/\${PROG}\n";
print MAKEFILE "\t\${", &LanguageCompiler($ARGV[1], @srcs);
print MAKEFILE "} \${LDFLAGS} $out\$@ \${OBJS} \${LIBS}\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \${OBJS} \${MODS}\n\n";
#
# make deinstall
#
print MAKEFILE "deinstall:\n";
print MAKEFILE "\trm -f \${OBJS} \${MODS} \${PROG}\n\n";
#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \${SUFFIXES} .f90 .F90 .c .obj\n\n";
#
# .f90 -> .o
#
print MAKEFILE ".f90.$o:\n";
print MAKEFILE "\t\${F90} \${F90FLAGS} \$<\n\n";
print MAKEFILE ".F90.$o:\n";
print MAKEFILE "\t\${F90} \${F90FLAGS} \$<\n\n";
print MAKEFILE ".c.$o:\n";
print MAKEFILE "\t\${CC} \${CFLAGS} \$<\n\n";
#
@fppFiles = qw();
foreach $fppFile (@fppFiles) {
  print MAKEFILE "$fppFile.obj: $fppFile.f90\n";
  print MAKEFILE "\t\${F90} \${F90FLAGS} /fpp \$<\n";
}
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
&MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "f90";
         }
      }
   else {
      CASE: {
         grep(/\.f90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   if ($platform eq "CYGWIN_NT-4.0") {
   $compiler = "LINK";
   }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.$o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90 *.F90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.[fF]90$/.$o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90 *.F90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.[fF]90$/.$o/;
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         #
         # Cray F90 compiler
         #
         if ($compiler eq "cray") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               push(@modules, "-p", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         #
         # ParaSoft F90 compiler
         #
         if ($compiler eq "parasoft") {
            print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c ";
            foreach $depend (@dependencies) {
               $depend =~ s/\.o$/.f90/;
               push(@modules, "-module", $depend);
               }
            push(@modules, $file);
            &PrintWords(30, 1, @modules);
            print MAKEFILE "\n";
            undef @modules;
            }
         }
      }
   }
