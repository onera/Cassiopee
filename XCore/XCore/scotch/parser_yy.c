/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 2

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         scotchyyparse
#define yylex           scotchyylex
#define yyerror         scotchyyerror
#define yydebug         scotchyydebug
#define yynerrs         scotchyynerrs


/* Copy the first part of user declarations.  */
#line 1 "parser_yy.y" /* yacc.c:339  */

/* Copyright 2004,2007,2008,2011,2014,2018,2019,2021,2023 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : parser_yy.y                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the syntactic parser     **/
/**                which processes strategy strings.       **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to   : 13 jun 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to   : 27 feb 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to   : 01 oct 1998     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to   : 11 jun 2004     **/
/**                # Version 5.1  : from : 30 oct 2007     **/
/**                                 to   : 24 jul 2011     **/
/**                # Version 6.0  : from : 30 sep 2014     **/
/**                                 to   : 27 apr 2018     **/
/**                # Version 7.0  : from : 02 mar 2018     **/
/**                                 to   : 20 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SCOTCH_PARSER_YY

#include "module.h"
#include "common.h"

#undef INTEGER                                    /* In case someone defined them */
#undef DOUBLE

/*
**  The type and structure definitions.
*/

/* Helpful definitions */

typedef void * yyscan_t;                          /* Should have been exported by Flex in y.tab.h */
typedef void * YY_BUFFER_STATE;                   /* The same; Flex and Bison design is just crap */

/*
**  The defines and includes (bis).
*/

#include "parser.h"
#include "parser_yy.h"
#include "parser_ly.h"
#include "parser_ll.h"

#ifdef SCOTCH_DEBUG_PARSER3
#define YYDEBUG                     1
#endif /* SCOTCH_DEBUG_PARSER3 */

/*
**  The static and global definitions.
**  See also at the end of this file.
*/

/*+ Method token conversion array. +*/

extern unsigned int         parsermethtokentab[];  /* Pre-definition */

/*+ Function prototypes that should have been exported by Flex. +*/

YY_BUFFER_STATE scotchyy_scan_string (const char *, yyscan_t);


#line 179 "parser_yy.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser_ly.h".  */
#ifndef YY_SCOTCHYY_PARSER_LY_H_INCLUDED
# define YY_SCOTCHYY_PARSER_LY_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int scotchyydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    METHODNAME = 258,
    PARAMNAME = 259,
    VALCASE = 260,
    VALDOUBLE = 261,
    VALINT = 262,
    VALSTRING = 263,
    VALSTRAT = 264,
    VALPARAM = 265,
    VALTEST = 266
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 122 "parser_yy.y" /* yacc.c:355  */

  char                      CASEVAL;              /* Case value          */
  StratTest *               TEST;                 /* Test type           */
  StratTestType             TESTOP;               /* Relational type     */
  double                    DOUBLE;               /* Double-precision    */
  INT                       INTEGER;              /* Integer             */
  char                      STRING[PARSERSTRINGLEN]; /* Character string */
  struct {
    const StratTab *        tabl;                 /* Current tables    */
    Strat *                 strat;                /* Current method    */
    const StratParamTab *   param;                /* Current parameter */
  } SAVE;                                         /* Parameter type    */
  Strat *                   STRAT;                /* Strategy tree     */

#line 246 "parser_yy.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif

/* Location type.  */
#if ! defined YYLTYPE && ! defined YYLTYPE_IS_DECLARED
typedef struct YYLTYPE YYLTYPE;
struct YYLTYPE
{
  int first_line;
  int first_column;
  int last_line;
  int last_column;
};
# define YYLTYPE_IS_DECLARED 1
# define YYLTYPE_IS_TRIVIAL 1
#endif



int scotchyyparse (void * scanptr, ParserEnv * penvptr);

#endif /* !YY_SCOTCHYY_PARSER_LY_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 276 "parser_yy.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL \
             && defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
  YYLTYPE yyls_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE) + sizeof (YYLTYPE)) \
      + 2 * YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  13
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   93

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  31
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  37
/* YYNRULES -- Number of rules.  */
#define YYNRULES  65
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  93

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   266

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    24,     2,     2,     2,    30,    23,     2,
      17,    18,    29,    27,    21,    28,     2,    13,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    16,    15,
      25,    22,    26,    14,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    19,    12,    20,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   162,   162,   168,   186,   189,   191,   206,   224,   228,
     232,   228,   256,   259,   264,   269,   273,   277,   276,   339,
     343,   339,   347,   350,   351,   355,   354,   393,   428,   443,
     458,   475,   475,   499,   508,   511,   529,   532,   550,   553,
     569,   573,   576,   595,   599,   603,   609,   625,   628,   632,
     638,   654,   657,   663,   679,   682,   688,   692,   693,   696,
     711,   728,   769,   773,   776,   780
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "METHODNAME", "PARAMNAME", "VALCASE",
  "VALDOUBLE", "VALINT", "VALSTRING", "VALSTRAT", "VALPARAM", "VALTEST",
  "'|'", "'/'", "'?'", "';'", "':'", "'('", "')'", "'{'", "'}'", "','",
  "'='", "'&'", "'!'", "'<'", "'>'", "'+'", "'-'", "'*'", "'%'", "$accept",
  "STRAT", "STRATSELECT", "STRATEMPTY", "STRATCONCAT", "STRATTEST", "$@1",
  "$@2", "STRATTESTELSE", "STRATGROUP", "STRATMETHOD", "$@3",
  "METHODPARAM", "$@4", "$@5", "PARAMLIST", "PARAMPARAM", "@6", "PARAMVAL",
  "@7", "TEST", "TESTOR", "TESTAND", "TESTNOT", "TESTREL", "TESTRELOP",
  "TESTEXPR1", "TESTEXPR1OP", "TESTEXPR2", "TESTEXPR2OP", "TESTEXPR3",
  "TESTEXPR3OP", "TESTEXPR4", "TESTVAL", "TESTVAR", "VALSDOUBLE",
  "VALSINT", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   124,    47,    63,    59,    58,    40,    41,   123,
     125,    44,    61,    38,    33,    60,    62,    43,    45,    42,
      37
};
# endif

#define YYPACT_NINF -33

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-33)))

#define YYTABLE_NINF -32

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int8 yypact[] =
{
       5,   -33,     5,    10,     7,   -33,     5,   -33,     4,   -33,
     -33,    19,    40,   -33,     5,   -33,    27,   -33,    20,   -33,
     -33,   -33,   -33,   -33,    27,    27,   -33,   -33,   -33,    28,
       3,   -33,   -33,    47,    -1,    33,    13,   -33,   -33,   -33,
     -33,   -33,    44,    41,    39,   -33,    42,    27,    27,   -33,
     -33,   -33,    43,    43,   -33,   -33,   -33,    43,   -33,    43,
     -33,    55,   -33,   -33,   -33,     5,     3,   -33,    43,    -4,
      33,    13,   -33,    46,    44,    57,    -9,    14,     8,   -33,
     -33,     5,    63,   -33,   -33,   -33,   -33,     5,   -33,   -33,
       7,   -33,     7
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       6,    17,     6,     0,     2,     4,     5,     8,     0,    12,
      16,    22,     0,     1,     6,     7,     0,    18,     0,    15,
       3,    61,    63,    65,     0,     0,    48,    49,    10,    34,
      36,    38,    41,     0,     0,    47,    51,    54,    57,    58,
      59,    60,     0,     0,     0,    39,     0,     0,     0,    44,
      43,    45,     0,     0,    62,    64,    52,     0,    55,     0,
      25,    20,    24,    40,    56,     6,    35,    37,     0,    42,
      46,    50,    53,     0,     0,     0,    14,     0,     0,    23,
      21,     6,     0,    33,    27,    30,    26,     6,    28,    29,
      13,    11,    32
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -33,   -33,    -2,    66,   -33,    75,   -33,   -33,   -33,   -33,
     -33,   -33,   -33,   -33,   -33,   -33,     9,   -33,   -33,   -33,
     -33,    58,    37,   -21,   -33,   -33,   -22,   -32,    34,   -33,
      29,   -33,    30,   -33,   -33,    12,    15
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     3,     4,     5,     6,     7,     8,    46,    82,     9,
      10,    11,    17,    18,    75,    61,    62,    73,    86,    87,
      28,    29,    30,    31,    32,    52,    33,    34,    35,    57,
      36,    59,    37,    38,    39,    40,    41
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      12,    53,    44,    14,    45,    54,    55,    81,     1,    83,
      13,   -31,    53,    84,    22,    23,    85,    16,    -9,    14,
     -31,   -31,     2,    26,    27,   -31,    48,    67,   -31,   -31,
      69,    21,    64,    22,    23,    26,    27,    53,   -19,    42,
      47,    26,    27,    58,    24,    53,    77,    21,    60,    22,
      23,    25,    14,    47,    26,    27,    65,    64,    19,    63,
      68,    49,    56,    76,    50,    51,    26,    27,    78,    49,
      26,    27,    50,    51,    26,    27,    74,    80,    91,    90,
      20,    15,    43,    79,    66,    92,    71,    70,     0,    72,
      88,     0,     0,    89
};

static const yytype_int8 yycheck[] =
{
       2,    33,    24,    12,    25,     6,     7,    16,     3,     1,
       0,     3,    44,     5,     6,     7,     8,    13,    13,    12,
      12,    13,    17,    27,    28,    17,    23,    48,    20,    21,
      52,     4,    18,     6,     7,    27,    28,    69,    19,    19,
      12,    27,    28,    30,    17,    77,    68,     4,     4,     6,
       7,    24,    12,    12,    27,    28,    14,    18,    18,    18,
      17,    22,    29,    65,    25,    26,    27,    28,    22,    22,
      27,    28,    25,    26,    27,    28,    21,    20,    15,    81,
      14,     6,    24,    74,    47,    87,    57,    53,    -1,    59,
      78,    -1,    -1,    78
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,    17,    32,    33,    34,    35,    36,    37,    40,
      41,    42,    33,     0,    12,    36,    13,    43,    44,    18,
      34,     4,     6,     7,    17,    24,    27,    28,    51,    52,
      53,    54,    55,    57,    58,    59,    61,    63,    64,    65,
      66,    67,    19,    52,    57,    54,    38,    12,    23,    22,
      25,    26,    56,    58,     6,     7,    29,    60,    30,    62,
       4,    46,    47,    18,    18,    14,    53,    54,    17,    57,
      59,    61,    63,    48,    21,    45,    33,    57,    22,    47,
      20,    16,    39,     1,     5,     8,    49,    50,    66,    67,
      33,    15,    33
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    31,    32,    33,    33,    34,    34,    35,    35,    37,
      38,    36,    36,    39,    39,    40,    40,    42,    41,    44,
      45,    43,    43,    46,    46,    48,    47,    49,    49,    49,
      49,    50,    49,    49,    51,    52,    52,    53,    53,    54,
      54,    54,    55,    56,    56,    56,    57,    57,    58,    58,
      59,    59,    60,    61,    61,    62,    63,    63,    63,    64,
      64,    65,    66,    66,    67,    67
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     3,     1,     1,     0,     2,     1,     0,
       0,     8,     1,     2,     0,     3,     1,     0,     3,     0,
       0,     5,     0,     3,     1,     0,     4,     1,     1,     1,
       1,     0,     2,     1,     1,     3,     1,     3,     1,     2,
       3,     1,     3,     1,     1,     1,     3,     1,     1,     1,
       3,     1,     1,     3,     1,     1,     3,     1,     1,     1,
       1,     1,     2,     1,     2,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (&yylloc, scanptr, penvptr, YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                                \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;        \
          (Current).first_column = YYRHSLOC (Rhs, 1).first_column;      \
          (Current).last_line    = YYRHSLOC (Rhs, N).last_line;         \
          (Current).last_column  = YYRHSLOC (Rhs, N).last_column;       \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).first_line   = (Current).last_line   =              \
            YYRHSLOC (Rhs, 0).last_line;                                \
          (Current).first_column = (Current).last_column =              \
            YYRHSLOC (Rhs, 0).last_column;                              \
        }                                                               \
    while (0)
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K])


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL

/* Print *YYLOCP on YYO.  Private, do not rely on its existence. */

YY_ATTRIBUTE_UNUSED
static unsigned
yy_location_print_ (FILE *yyo, YYLTYPE const * const yylocp)
{
  unsigned res = 0;
  int end_col = 0 != yylocp->last_column ? yylocp->last_column - 1 : 0;
  if (0 <= yylocp->first_line)
    {
      res += YYFPRINTF (yyo, "%d", yylocp->first_line);
      if (0 <= yylocp->first_column)
        res += YYFPRINTF (yyo, ".%d", yylocp->first_column);
    }
  if (0 <= yylocp->last_line)
    {
      if (yylocp->first_line < yylocp->last_line)
        {
          res += YYFPRINTF (yyo, "-%d", yylocp->last_line);
          if (0 <= end_col)
            res += YYFPRINTF (yyo, ".%d", end_col);
        }
      else if (0 <= end_col && yylocp->first_column < end_col)
        res += YYFPRINTF (yyo, "-%d", end_col);
    }
  return res;
 }

#  define YY_LOCATION_PRINT(File, Loc)          \
  yy_location_print_ (File, &(Loc))

# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value, Location, scanptr, penvptr); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, void * scanptr, ParserEnv * penvptr)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (yylocationp);
  YYUSE (scanptr);
  YYUSE (penvptr);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, YYLTYPE const * const yylocationp, void * scanptr, ParserEnv * penvptr)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  YY_LOCATION_PRINT (yyoutput, *yylocationp);
  YYFPRINTF (yyoutput, ": ");
  yy_symbol_value_print (yyoutput, yytype, yyvaluep, yylocationp, scanptr, penvptr);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, YYLTYPE *yylsp, int yyrule, void * scanptr, ParserEnv * penvptr)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                       , &(yylsp[(yyi + 1) - (yynrhs)])                       , scanptr, penvptr);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, yylsp, Rule, scanptr, penvptr); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, YYLTYPE *yylocationp, void * scanptr, ParserEnv * penvptr)
{
  YYUSE (yyvaluep);
  YYUSE (yylocationp);
  YYUSE (scanptr);
  YYUSE (penvptr);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void * scanptr, ParserEnv * penvptr)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

/* Location data for the lookahead symbol.  */
static YYLTYPE yyloc_default
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
  = { 1, 1, 1, 1 }
# endif
;
YYLTYPE yylloc = yyloc_default;

    /* Number of syntax errors so far.  */
    int yynerrs;

    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.
       'yyls': related to locations.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls;
    YYLTYPE *yylsp;

    /* The locations where the error started and ended.  */
    YYLTYPE yyerror_range[3];

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;
  YYLTYPE yyloc;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N), yylsp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yylsp = yyls = yylsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

/* User initialization code.  */
#line 113 "parser_yy.y" /* yacc.c:1429  */
{
  yylloc.cobenum =
  yylloc.libenum =
  yylloc.coennum =
  yylloc.liennum = 1;
  yylloc.tebeptr =
  yylloc.teenptr = penvptr->textptr;
}

#line 1352 "parser_yy.c" /* yacc.c:1429  */
  yylsp[0] = yylloc;
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;
        YYLTYPE *yyls1 = yyls;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yyls1, yysize * sizeof (*yylsp),
                    &yystacksize);

        yyls = yyls1;
        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
        YYSTACK_RELOCATE (yyls_alloc, yyls);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;
      yylsp = yyls + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex (&yylval, &yylloc, scanptr);
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END
  *++yylsp = yylloc;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];

  /* Default location.  */
  YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 163 "parser_yy.y" /* yacc.c:1646  */
    {
                penvptr->straptr = ((yyvsp[0].STRAT));          /* Save pointer to root of tree */
              }
#line 1543 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 3:
#line 169 "parser_yy.y" /* yacc.c:1646  */
    {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (1)");
                  stratExit  ((yyvsp[-2].STRAT));
                  stratExit  ((yyvsp[0].STRAT));
                  YYABORT;
                }

                straptr->tabl                 = penvptr->stratab;
                straptr->type                 = STRATNODESELECT;
                straptr->data.select.strat[0] = ((yyvsp[-2].STRAT));
                straptr->data.select.strat[1] = ((yyvsp[0].STRAT));

                ((yyval.STRAT)) = straptr;
              }
#line 1565 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 6:
#line 191 "parser_yy.y" /* yacc.c:1646  */
    {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (2)");
                  YYABORT;
                }

                straptr->tabl = penvptr->stratab;
                straptr->type = STRATNODEEMPTY;

                ((yyval.STRAT)) = straptr;
              }
#line 1583 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 7:
#line 207 "parser_yy.y" /* yacc.c:1646  */
    {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (3)");
                  stratExit  ((yyvsp[-1].STRAT));
                  stratExit  ((yyvsp[0].STRAT));
                  YYABORT;
                }

                straptr->tabl                 = penvptr->stratab;
                straptr->type                 = STRATNODECONCAT;
                straptr->data.concat.strat[0] = ((yyvsp[-1].STRAT));
                straptr->data.concat.strat[1] = ((yyvsp[0].STRAT));

                ((yyval.STRAT)) = straptr;
              }
#line 1605 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 9:
#line 228 "parser_yy.y" /* yacc.c:1646  */
    {
                PARSERLLBEGIN (VALTEST);          /* Parse parameter tokens */
              }
#line 1613 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 10:
#line 232 "parser_yy.y" /* yacc.c:1646  */
    {
                PARSERLLBEGIN (VALSTRAT);         /* Parse strategy tokens */
              }
#line 1621 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 11:
#line 236 "parser_yy.y" /* yacc.c:1646  */
    {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (4)");
                  stratExit  ((yyvsp[-2].STRAT));
                  if (((yyvsp[-1].STRAT)) != NULL)
                    stratExit ((yyvsp[-1].STRAT));
                  stratTestExit ((yyvsp[-5].TEST));
                  YYABORT;
                }

                straptr->tabl               = penvptr->stratab;
                straptr->type               = STRATNODECOND;
                straptr->data.cond.test     = ((yyvsp[-5].TEST));
                straptr->data.cond.strat[0] = ((yyvsp[-2].STRAT));
                straptr->data.cond.strat[1] = ((yyvsp[-1].STRAT));

                ((yyval.STRAT)) = straptr;
              }
#line 1646 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 13:
#line 260 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.STRAT)) = ((yyvsp[0].STRAT));
              }
#line 1654 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 14:
#line 264 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.STRAT)) = NULL;
              }
#line 1662 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 15:
#line 270 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.STRAT)) = ((yyvsp[-1].STRAT));
              }
#line 1670 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 17:
#line 277 "parser_yy.y" /* yacc.c:1646  */
    {
                Strat *             straptr;
                int                 methnum;
                int                 methlen;
                StratMethodTab *    methtab;
                int                 i, j;

                methnum =
                methlen = 0;                      /* No method recognized yet     */
                methtab = penvptr->stratab->methtab; /* Point to the method table */
                for (i = 0; methtab[i].name != NULL; i ++) {
                  if ((strncasecmp (((yyvsp[0].STRING)),         /* Find longest matching code name */
                       methtab[i].name,
                       j = strlen (methtab[i].name)) == 0) &&
                      (j > methlen)) {
                    methnum = methtab[i].meth;
                    methlen = j;
                  }
                }
                if (methlen == 0) {               /* If method name not known */
                  errorPrint ("stratParserParse: invalid method name \"%s\", line %d, column %d, at \"%s\"",
                              ((yyvsp[0].STRING)), ((yylsp[0])).libenum, ((yylsp[0])).cobenum, ((yylsp[0])).tebeptr);
                  YYABORT;
                }
                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (5)");
                  YYABORT;
                }

                straptr->tabl             = penvptr->stratab;
                straptr->type             = STRATNODEMETHOD;
                straptr->data.method.meth = methnum; /* Set method type        */
                if (methtab[methnum].data != NULL) /* If default values exist  */
                  memcpy (&straptr->data.method.data, /* Set values to default */
                          methtab[methnum].data,
                          sizeof (StratNodeMethodData));

                penvptr->straptr = straptr;       /* Structure available for parameter processing */
              }
#line 1714 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 18:
#line 317 "parser_yy.y" /* yacc.c:1646  */
    {
                StratParamTab *     paratab;
                int                 paraidx;

                paratab = penvptr->stratab->paratab; /* Point to the parameter table */
                for (paraidx = 0; paratab[paraidx].name != NULL; paraidx ++) {
                  if ((paratab[paraidx].meth == penvptr->straptr->data.method.meth) && /* If a strategy parameter found for this method */
                      (paratab[paraidx].type == STRATPARAMSTRAT)) {
                    if (*((Strat **) ((byte *) &penvptr->straptr->data.method.data + /* And this parameter has not been set */
                        (paratab[paraidx].dataofft - paratab[paraidx].database))) == NULL)
                      errorPrintW ("stratParserParse: strategy parameter \"%s\" of method \"%s\" not set, line %d, column %d, before \"%s\"",
                                   paratab[paraidx].name, penvptr->stratab->methtab[penvptr->straptr->data.method.meth].name,
                                   ((yylsp[-1])).libenum, ((yylsp[-1])).cobenum, ((yylsp[-1])).tebeptr);
                  }
                }

                ((yyval.STRAT)) = penvptr->straptr;          /* Return current structure */
                penvptr->straptr = NULL;          /* No current structure     */
              }
#line 1738 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 19:
#line 339 "parser_yy.y" /* yacc.c:1646  */
    {
                PARSERLLBEGIN (VALPARAM);         /* Parse parameter tokens */
              }
#line 1746 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 20:
#line 343 "parser_yy.y" /* yacc.c:1646  */
    {
                PARSERLLBEGIN (VALSTRAT);         /* Parse strategy tokens */
              }
#line 1754 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 25:
#line 355 "parser_yy.y" /* yacc.c:1646  */
    {
                int               paraidx;
                int               paralen;
                StratParamTab *   paratab;
                int               i, j;

                paraidx =
                paralen = 0;                      /* No parameter recognized yet     */
                paratab = penvptr->stratab->paratab; /* Point to the parameter table */
                for (i = 0; paratab[i].name != NULL; i ++) {
                  if ((paratab[i].meth == penvptr->straptr->data.method.meth) &&
                      (strncasecmp (((yyvsp[0].STRING)),         /* Find longest matching parameter name */
                                    paratab[i].name,
                                    j = strlen (paratab[i].name)) == 0) &&
                      (j > paralen)) {
                    paraidx = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid method parameter name \"%s\", line %d, column %d, before \"%s\"",
                              ((yyvsp[0].STRING)), ((yylsp[0])).libenum, ((yylsp[0])).cobenum, ((yylsp[0])).tebeptr);
                  YYABORT;
                }

                ((yyval.SAVE)).tabl = penvptr->stratab; /* Save current strategy tables   */
                penvptr->paraptr = &paratab[paraidx]; /* Save current parameter value */
                PARSERLLBEGIN (stratmethtokentab[penvptr->paraptr->type & ~STRATPARAMDEPRECATED]); /* Get non-deprecated type */
                if (penvptr->paraptr->type == STRATPARAMSTRAT) /* If parameter is a strategy                                  */
                  penvptr->stratab = (StratTab *) penvptr->paraptr->datasltr; /* Use new strategy tables                      */
              }
#line 1790 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 26:
#line 387 "parser_yy.y" /* yacc.c:1646  */
    {
                PARSERLLBEGIN (VALPARAM);         /* Go-on reading parameters          */
                penvptr->stratab = ((yyvsp[-2].SAVE)).tabl; /* Restore current strategy tables */
              }
#line 1799 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 27:
#line 394 "parser_yy.y" /* yacc.c:1646  */
    {
                char              c;              /* Character read             */
                char *            p;              /* Pointer to selector string */
                int               i;              /* Index in selector string   */

                if ((penvptr->paraptr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
                  c = ((yyvsp[0].CASEVAL));                       /* First, use char as is */
                  for (p = (char *) penvptr->paraptr->datasltr, i = 0;
                       (*p != '\0') && (*p != c);
                       p ++, i ++) ;
                  if (*p == '\0') {               /* Char was not found         */
                    c = tolower (c);              /* Convert char to lower case */
                    for (p = (char *) penvptr->paraptr->datasltr, i = 0;
                         (*p != '\0') && (*p != c);
                         p ++, i ++) ;
                    if (*p == '\0') {
                      errorPrint ("stratParserParse: invalid method parameter switch \"%s=%c\", line %d, column %d, before \"%s\"",
                                  penvptr->paraptr->name, ((yyvsp[0].CASEVAL)), ((yylsp[0])).libenum, ((yylsp[0])).cobenum, ((yylsp[0])).tebeptr);
                      YYABORT;
                    }
                  }

#ifdef SCOTCH_DEBUG_PARSER2
                  if ((penvptr->paraptr->dataofft - penvptr->paraptr->database + sizeof (int)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (1)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((int *) ((byte *) &penvptr->straptr->data.method.data +
                             (penvptr->paraptr->dataofft -
                              penvptr->paraptr->database))) = i;
                }
              }
#line 1838 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 28:
#line 429 "parser_yy.y" /* yacc.c:1646  */
    {
                if ((penvptr->paraptr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((penvptr->paraptr->dataofft - penvptr->paraptr->database + sizeof (double)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (2)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((double *) ((byte *) penvptr->straptr->data.method.data +
                                (penvptr->paraptr->dataofft -
                                 penvptr->paraptr->database))) = ((yyvsp[0].DOUBLE));
                }
              }
#line 1857 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 29:
#line 444 "parser_yy.y" /* yacc.c:1646  */
    {
                if ((penvptr->paraptr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((penvptr->paraptr->dataofft - penvptr->paraptr->database + sizeof (INT)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (3)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((INT *) ((byte *) &penvptr->straptr->data.method.data +
                             (penvptr->paraptr->dataofft -
                              penvptr->paraptr->database))) = (INT) ((yyvsp[0].INTEGER));
                }
              }
#line 1876 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 30:
#line 459 "parser_yy.y" /* yacc.c:1646  */
    {
                if ((penvptr->paraptr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((penvptr->paraptr->dataofft - penvptr->paraptr->database + strlen ((yyvsp[0].STRING)) + 1) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (4)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  strcpy ((char *) ((byte *) &penvptr->straptr->data.method.data +
                                    (penvptr->paraptr->dataofft -
                                     penvptr->paraptr->database)),
                          ((yyvsp[0].STRING)));
                }
              }
#line 1896 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 31:
#line 475 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.SAVE)).strat = penvptr->straptr;
                ((yyval.SAVE)).param = penvptr->paraptr;
                penvptr->straptr = NULL;
                penvptr->paraptr = NULL;
              }
#line 1907 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 32:
#line 482 "parser_yy.y" /* yacc.c:1646  */
    {
                penvptr->straptr = ((yyvsp[-1].SAVE)).strat; /* Restore current method    */
                penvptr->paraptr = ((yyvsp[-1].SAVE)).param; /* Restore current parameter */

                if ((penvptr->paraptr->type & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if ((penvptr->paraptr->dataofft - penvptr->paraptr->database + sizeof (Strat *)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (5)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((Strat **) ((byte *) &penvptr->straptr->data.method.data +
                                (penvptr->paraptr->dataofft -
                                 penvptr->paraptr->database))) = ((yyvsp[0].STRAT));
                }
              }
#line 1929 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 33:
#line 500 "parser_yy.y" /* yacc.c:1646  */
    {
                errorPrint ("stratParserParse: invalid value for parameter \"%s\" of method \"%s\", line %d, column %d, before \"%s\"",
                            penvptr->paraptr->name, penvptr->straptr->tabl->methtab[penvptr->straptr->data.method.meth].name,
                            ((yylsp[0])).libenum, ((yylsp[0])).cobenum, ((yylsp[0])).tebeptr);
                YYABORT;
              }
#line 1940 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 35:
#line 512 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (6)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }

                testptr->typetest     = STRATTESTOR;
                testptr->typenode     = STRATPARAMLOG;
                testptr->data.test[0] = ((yyvsp[-2].TEST));
                testptr->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 1962 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 37:
#line 533 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (7)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }

                testptr->typetest     = STRATTESTAND;
                testptr->typenode     = STRATPARAMLOG;
                testptr->data.test[0] = ((yyvsp[-2].TEST));
                testptr->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 1984 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 39:
#line 554 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (8)");
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }

                testptr->typetest     = STRATTESTNOT;
                testptr->typenode     = STRATPARAMLOG;
                testptr->data.test[0] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 2004 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 40:
#line 570 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TEST)) = ((yyvsp[-1].TEST));
              }
#line 2012 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 42:
#line 577 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (9)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }
                testptr->typetest     = ((yyvsp[-1].TESTOP));
                testptr->typenode     = STRATPARAMLOG;
                testptr->data.test[0] = ((yyvsp[-2].TEST));
                testptr->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 2033 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 43:
#line 596 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTLT;
              }
#line 2041 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 44:
#line 600 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTEQ;
              }
#line 2049 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 45:
#line 604 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTGT;
              }
#line 2057 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 46:
#line 610 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (10)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }
                testptr->typetest     = ((yyvsp[-1].TESTOP));
                testptr->data.test[0] = ((yyvsp[-2].TEST));
                testptr->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 2077 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 48:
#line 629 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTADD;
              }
#line 2085 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 49:
#line 633 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTSUB;
              }
#line 2093 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 50:
#line 639 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  errorPrint    ("stratParserParse: out of memory (11)");
                  YYABORT;
                }
                testptr->typetest     = ((yyvsp[-1].TESTOP));
                testptr->data.test[0] = ((yyvsp[-2].TEST));
                testptr->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 2113 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 52:
#line 658 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTMUL;
              }
#line 2121 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 53:
#line 664 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (12)");
                  stratTestExit ((yyvsp[-2].TEST));
                  stratTestExit ((yyvsp[0].TEST));
                  YYABORT;
                }
                testptr->typetest     = ((yyvsp[-1].TESTOP));
                testptr->data.test[0] = ((yyvsp[-2].TEST));
                testptr->data.test[1] = ((yyvsp[0].TEST));

                ((yyval.TEST)) = testptr;
              }
#line 2141 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 55:
#line 683 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TESTOP)) = STRATTESTMOD;
              }
#line 2149 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 56:
#line 689 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.TEST)) = ((yyvsp[-1].TEST));
              }
#line 2157 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 59:
#line 697 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (13)");
                  YYABORT;
                }

                testptr->typetest        = STRATTESTVAL;
                testptr->typenode        = STRATPARAMDOUBLE;
                testptr->data.val.valdbl = ((yyvsp[0].DOUBLE));

                ((yyval.TEST)) = testptr;
              }
#line 2176 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 60:
#line 712 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (14)");
                  YYABORT;
                }

                testptr->typetest        = STRATTESTVAL;
                testptr->typenode        = STRATPARAMINT;
                testptr->data.val.valint = ((yyvsp[0].INTEGER));

                ((yyval.TEST)) = testptr;
              }
#line 2195 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 61:
#line 729 "parser_yy.y" /* yacc.c:1646  */
    {
                StratTest *       testptr;
                StratParamTab *   condtab;
                int               para;
                int               paralen;
                int               i, j;

                para    =
                paralen = 0;                      /* No parameter recognized yet */
                condtab = penvptr->stratab->condtab; /* Point to parameter table */
                for (i = 0; condtab[i].name != NULL; i ++) {
                  if ((strncasecmp (((yyvsp[0].STRING)),         /* Find longest matching parameter name */
                                    condtab[i].name,
                                    j = strlen (condtab[i].name)) == 0) &&
                      (j > paralen)) {
                    para    = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid graph parameter name \"%s\", line %d, column %d, before \"%s\"",
                              ((yyvsp[0].STRING)), ((yylsp[0])).libenum, ((yylsp[0])).cobenum, ((yylsp[0])).tebeptr);
                  YYABORT;
                }

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (15)");
                  YYABORT;
                }

                testptr->typetest          = STRATTESTVAR;
                testptr->typenode          = condtab[para].type;
                testptr->data.var.datatab  = penvptr->stratab;
                testptr->data.var.datadisp = condtab[para].dataofft -
                                             condtab[para].database;

                ((yyval.TEST)) = testptr;
              }
#line 2238 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 62:
#line 770 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.DOUBLE)) = (((yyvsp[-1].TESTOP)) == STRATTESTSUB) ? - ((yyvsp[0].DOUBLE)) : ((yyvsp[0].DOUBLE));
              }
#line 2246 "parser_yy.c" /* yacc.c:1646  */
    break;

  case 64:
#line 777 "parser_yy.y" /* yacc.c:1646  */
    {
                ((yyval.INTEGER)) = (((yyvsp[-1].TESTOP)) == STRATTESTSUB) ? - ((yyvsp[0].INTEGER)) : ((yyvsp[0].INTEGER));
              }
#line 2254 "parser_yy.c" /* yacc.c:1646  */
    break;


#line 2258 "parser_yy.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;
  *++yylsp = yyloc;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (&yylloc, scanptr, penvptr, YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (&yylloc, scanptr, penvptr, yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }

  yyerror_range[1] = yylloc;

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval, &yylloc, scanptr, penvptr);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  yyerror_range[1] = yylsp[1-yylen];
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      yyerror_range[1] = *yylsp;
      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp, yylsp, scanptr, penvptr);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  yyerror_range[2] = yylloc;
  /* Using YYLLOC is tempting, but would change the location of
     the lookahead.  YYLOC is available though.  */
  YYLLOC_DEFAULT (yyloc, yyerror_range, 2);
  *++yylsp = yyloc;

  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (&yylloc, scanptr, penvptr, YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval, &yylloc, scanptr, penvptr);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, yylsp, scanptr, penvptr);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 783 "parser_yy.y" /* yacc.c:1906  */


/*
**  The static and global definitions (bis).
**  These are put at the end of the file because
**  the token values that they use are not yet
**  defined in the first section of the file.
*/

unsigned int                stratmethtokentab[] = { /* Table for parameter/token type conversion */
                              VALCASE,
                              VALDOUBLE,
                              VALINT,
                              -1,                 /* No logical parameters */
                              VALSTRAT,
                              VALSTRING,
                              -1                  /* One more value to detect array overflow */
                            };

/************************************/
/*                                  */
/* These routines drive the parser. */
/*                                  */
/************************************/

/* This routine is the entry point for
** the strategy parser.
** It returns:
** - !NULL  : pointer to the strategy.
** - NULL   : on error.
*/

Strat *
stratParserParse (
const StratTab * const      stratab,              /*+ Pointer to parsing tables +*/
const char * const          textptr)              /*+ Strategy string to parse  +*/
{
  YY_BUFFER_STATE     buffdat;
  ParserEnv           penvdat;                    /* Parser environment    */
  yyscan_t            scandat;                    /* Pointer to lex memory */
  int                 o;

  penvdat.stratab = stratab;                      /* Point to the parsing tables             */
  penvdat.straptr = NULL;                         /* Clear up the temporary strategy pointer */
  penvdat.textptr = textptr;                      /* Initialize the lexical parser           */

  if (scotchyylex_init (&scandat) != 0) {
    errorPrint ("stratParserParse: cannot initialize reentrant parser");
    return     (NULL);
  }
  buffdat = scotchyy_scan_string (textptr, scandat); /* Let's hope nothing breaks; error management in flex is just crap */
  scotchyy_switch_to_buffer (buffdat, scandat);

  o = yyparse (scandat, &penvdat);                /* Parse the strategy string */

  scotchyy_delete_buffer (buffdat, scandat);
  scotchyylex_destroy (scandat);

  if (o != 0) {
    if (penvdat.straptr != NULL)
      stratExit (penvdat.straptr);
    return (NULL);
  }

  return (penvdat.straptr);                       /* Return strategy pointer */
}

/* This routine displays the parser error message.
** It returns:
** - void  : in all cases.
*/

static
void
yyerror (
const ParserLocation * const  plocptr,            /*+ Scan location +*/
void * const                  scanptr,            /*+ Not used      +*/
const ParserEnv * const       penvptr,            /*+ Not used      +*/
const char * const            mesgptr)            /*+ Not used      +*/
{
  errorPrint ("stratParserParse: invalid strategy string, line %d, column %d, at \"%s\"",
              plocptr->libenum, plocptr->cobenum, plocptr->tebeptr);
}
