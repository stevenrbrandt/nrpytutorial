# An attempt at a WaveToy using the NRPy+ infrastructure

import os, re
from datetime import date
from sympy import symbols, Function, diff
import grid
import NRPy_param_funcs as par
import finite_difference as fin
from outputC import indent_Ccode, add_to_Cfunction_dict, outCfunction, construct_NRPy_function_prototypes_h, construct_NRPy_Cfunctions
import sympy
from outputC import lhrh
import grid as gri
import NRPy_param_funcs as par
from safewrite import SafeWrite

today = date.today().strftime("%B %d, %Y")

makefile_init="""
# Main make.code.defn file for thorn Foo

# Subdirectories containing source files
SUBDIRS =

# Source files in this directory
SRCS ="""

doc_init=r"""
% *======================================================================*
%  Cactus Thorn template for ThornGuide documentation
%  Author: {author}
%  Date: {date}
%  $Header$
%
%  Thorn documentation in the latex file doc/documentation.tex
%  will be included in ThornGuides built with the Cactus make system.
%  The scripts employed by the make system automatically include
%  pages about variables, parameters and scheduling parsed from the
%  relevant thorn CCL files.
%
%  This template contains guidelines which help to assure that your
%  documentation will be correctly added to ThornGuides. More
%  information is available in the Cactus UsersGuide.
%
%  Guidelines:
%   - Do not change anything before the line
%       % START CACTUS THORNGUIDE",
%     except for filling in the title, author, date, etc. fields.
%        - Each of these fields should only be on ONE line.
%        - Author names should be separated with a \\ or a comma.
%   - You can define your own macros, but they must appear after
%     the START CACTUS THORNGUIDE line, and must not redefine standard
%     latex commands.
%   - To avoid name clashes with other thorns, 'labels', 'citations',
%     'references', and 'image' names should conform to the following
%     convention:
%       ARRANGEMENT_THORN_LABEL
%     For example, an image wave.eps in the arrangement CactusWave and
%     thorn WaveToyC should be renamed to CactusWave_WaveToyC_wave.eps
%   - Graphics should only be included using the graphicx package.
%     More specifically, with the "\includegraphics" command.  Do
%     not specify any graphic file extensions in your .tex file. This
%     will allow us to create a PDF version of the ThornGuide
%     via pdflatex.
%   - References should be included with the latex "\bibitem" command.
%   - Use \begin{abstract}...\end{abstract} instead of \abstract{...}
%   - Do not use \appendix, instead include any appendices you need as
%     standard sections.
%   - For the benefit of our Perl scripts, and for future extensions,
%     please use simple latex.
%
% *======================================================================*
%
% Example of including a graphic image:
%    \begin{figure}[ht]
% 	\begin{center}
%    	   \includegraphics[width=6cm]{MyArrangement_MyThorn_MyFigure}
% 	\end{center}
% 	\caption{Illustration of this and that}
% 	\label{MyArrangement_MyThorn_MyLabel}
%    \end{figure}
%
% Example of using a label:
%   \label{MyArrangement_MyThorn_MyLabel}
%
% Example of a citation:
%    \cite{MyArrangement_MyThorn_Author99}
%
% Example of including a reference
%   \bibitem{MyArrangement_MyThorn_Author99}
%   {J. Author, {\em The Title of the Book, Journal, or periodical}, 1 (1999),
%   1--16. {\tt http://www.nowhere.com/}}
%
% *======================================================================*

% If you are using CVS use this line to give version information
% $Header$

\documentclass{article}

% Use the Cactus ThornGuide style file
% (Automatically used from Cactus distribution, if you have a
%  thorn without the Cactus Flesh download this from the Cactus
%  homepage at www.cactuscode.org)
\usepackage{../../../../doc/latex/cactus}

\begin{document}

% The author of the documentation
\author{{author} \textless sbrandt@cct.lsu.edu\textgreater}

% The title of the document (not necessarily the name of the Thorn)
\title{{thornname}}

% the date your document was last changed, if your document is in CVS,
% please use:
%    \date{$ $Date$ $}
% when using git instead record the commit ID:
%    \date{\gitrevision{<path-to-your-.git-directory>}}
\date{{date}}

\maketitle

% Do not delete next line
% START CACTUS THORNGUIDE

% Add all definitions used in this documentation here
%   \def\mydef etc

% Add an abstract for this thorn's documentation
\begin{abstract}

\end{abstract}

% The following sections are suggestive only.
% Remove them or add your own.

\section{Introduction}

\section{Physical System}

\section{Numerical Implementation}

\section{Using This Thorn}

\subsection{Obtaining This Thorn}

\subsection{Basic Usage}

\subsection{Special Behaviour}

\subsection{Interaction With Other Thorns}

\subsection{Examples}

\subsection{Support and Feedback}

\section{History}

\subsection{Thorn Source Code}

\subsection{Thorn Documentation}

\subsection{Acknowledgements}


\begin{thebibliography}{9}

\end{thebibliography}

% Do not delete next line
% END CACTUS THORNGUIDE

\end{document}
""".lstrip()

def typeof(*args):
    t = None
    for a in args:
        if a is None:
            continue
        ta = type(a)
        assert ta in [int, float, bool]
        if t is None:
            t = ta
        if t == int and ta == float:
            t = float
        elif t == float and ta == int:
            pass
        else:
            assert t == ta, f"{t} and {ta}"
    assert t is not None
    return t

class CactusFunc:
    def __init__(self, name, body, schedule_bin, doc):
        self.name = name
        self.body = body
        self.writegfs = set()
        self.readgfs = set()
        self.schedule_bin = schedule_bin
        self.doc = doc

class CactusSrc:
    def __init__(self, name):
        self.name = name
        self.funcs = []
        self.where = None
    def add_func(self, func, where):
        self.funcs += [func]
        self.where = where

class CactusThorn:

    def _add_src(self, src):
        if type(src) == str:
            self.src_files += [CactusSrc(src)]
        elif type(src) == CactusSrc:
            self.src_files += [src]
        else:
            assert False, "Bad type for src"

    def add_func(self, name, body, schedule_bin, doc, where='interior'):
        centering = None
        self._add_src(name + ".cc")
        writegfs = set()
        readgfs = set()
        has_fd = False
        if type(body)==list:
            new_body = []
            for item in body:
                assert type(item) == lhrh, "Equations should be stored in outputC.lhrh objects."
                writegfs.add(str(item.lhs))
                if hasattr(item.rhs, "free_symbols"):
                    for sym in item.rhs.free_symbols:
                        rdsym = str(sym)
                        # Check if the symbol name is a derivative
                        g = re.match(r'(.*)_dD+\d+$', rdsym)
                        if g:
                            readgfs.add(g.group(1))
                        else:
                            readgfs.add(rdsym)
                if type(item.lhs) == sympy.core.symbol.Symbol:
                    new_lhs=gri.gfaccess(varname=str(item.lhs))
                    new_body += [lhrh(lhs=new_lhs, rhs=item.rhs)]
                else:
                    new_body += [item]
            body = new_body
            kernel = fin.FD_outputC("returnstring",body)
            decl = ""
            if gri.ET_driver == "Carpet":
                body = f"""
                {decl}
                for(int i2=cctk_nghostzones[2];i2<cctk_lsh[2]-cctk_nghostzones[2];i2++) {{
                for(int i1=cctk_nghostzones[1];i1<cctk_lsh[1]-cctk_nghostzones[1];i1++) {{
                for(int i0=cctk_nghostzones[0];i0<cctk_lsh[0]-cctk_nghostzones[0];i0++) {{
                {kernel}
                }} }} }}
                """.strip()
                assert where == "interior", "Only interior suppprted right now"
            elif gri.ET_driver == "CarpetX":
                if where == 'everywhere':
                    wtag = 'all'
                elif where == 'interior':
                    wtag = 'int'
                elif where == 'boundary':
                    wtag = 'bnd'
                else:
                    assert False, "where should be in ['interior', 'everywhere', 'boundary']"
                #TODO: need to decide how to do loops with centering other than VVV
                body = f"""
                {decl}
                grid.loop_{wtag}_device<VVV_centered[0], VVV_centered[1], VVV_centered[2]>(
                grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {{
                {kernel}
                }});
                """.strip()
        elif type(body)==str:
            # Pass the body through literally
            pass
        else:
            assert False, "Body must be list or str"
        func = CactusFunc(name, body, schedule_bin, doc)
        self.src_files[-1].add_func(func, where)
        func.readgfs = readgfs
        func.writegfs = writegfs

    def declare_param(self, name, default, doc, vmin=None, vmax=None, options=None):
        self.params += [(name, default, doc, vmin, vmax, options)]
        ty = type(default)
        if ty == int:
            c_type = "int"
        elif ty == float:
            c_type = "REAL"
        else:
            raise Exception("Unkown parameter type: "+str(ty))
        par.Cparameters(c_type,self.thornname,[name],default)
        return symbols(name)

    def use_coords(self):
        return self.coords

    def register_gridfunctions(self, gtype, gf_names, centering=None, external_module=None):
        if external_module is not None:
            assert gtype == "EXTERNAL"
            thorn = external_module
        else:
            thorn = self.thornname
        for gf_name in gf_names:
            self.gf_names[gf_name] = thorn
            if centering is not None:
                self.centering[gf_name] = centering
            else:
                self.centering[gf_name] = "VVV"
        if external_module is not None:
            self.external_modules += [external_module]
        return grid.register_gridfunctions(gtype, gf_names, external_module=external_module)

    def __init__(self, arrangement, thornname, author=None, email=None, license='BSD'):
        self.ET_driver = grid.ET_driver
        self.gf_names = {}
        self.arrangement = arrangement
        self.thornname = thornname
        self.thorn_dir = os.path.join("arrangements", self.arrangement, self.thornname)

        self.src_files = []
        self.params = []
        self.centering = {}
        self.external_modules = []

        self.param_ccl = os.path.join(self.thorn_dir, "param.ccl")
        self.interface_ccl = os.path.join(self.thorn_dir, "interface.ccl")
        self.schedule_ccl = os.path.join(self.thorn_dir, "schedule.ccl")
        self.configuration_ccl = os.path.join(self.thorn_dir, "configuration.ccl")

        self.src_dir = os.path.join(self.thorn_dir, "src")
        self.makefile = os.path.join(self.src_dir, "make.code.defn")

        self.test_dir = os.path.join(self.thorn_dir, "test")
        self.test_ccl = os.path.join(self.test_dir, "test.ccl")

        self.doc_dir = os.path.join(self.thorn_dir, "doc")
        self.doc_tex = os.path.join(self.doc_dir, "documentation.tex")
        self.author = author
        self.email = email
        self.license = license

        par.set_parval_from_str("grid::DIM",3)
        #par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
        par.set_parval_from_str("grid::GridFuncMemAccess","ETK")

        if self.author is None or self.email is None:
            # Does gitconfig exist? Look up author and email there...
            home = os.environ["HOME"]
            git_config = os.path.join(home, ".gitconfig")
            if os.path.exists(git_config):
                import configparser
                config = configparser.ConfigParser(strict=False)
                config.read(git_config)
                if 'user' in config:
                    user = config['user']
                    if 'name' in user and self.author is None:
                        self.author = user['name']
                    if 'email' in user and self.email is None:
                        self.email = user['email']
        if self.author is None:
            self.author = os.environ['USER']
        if self.email is None:
            self.email = os.environ['USER']

    def get_full_name(self,gf_name):
        gfthorn = self.gf_names[gf_name]
        if gfthorn == self.thornname:
            return f"{gfthorn}::{gf_name}GF"
        else:
            return f"{gfthorn}::{gf_name}"

    def generate(self,dirname=None,cactus_config="sim",cactus_thornlist=None):
        assert self.ET_driver == grid.ET_driver
        cwd = None
        try:
            if dirname is not None:
                cwd = os.getcwd()
                os.chdir(dirname)
            assert os.path.isdir("arrangements"), "Please run this script from a Cactus root directory"
            os.makedirs(self.src_dir, exist_ok=True)
            os.makedirs(self.test_dir, exist_ok=True)
            os.makedirs(self.doc_dir, exist_ok=True)
            with SafeWrite(self.configuration_ccl) as fd:
                print(f"# Configuration definitions for thorn {self.thornname}",file=fd)
                if gri.ET_driver == "CarpetX":
                    print(f"REQUIRES CarpetX",file=fd)
            with SafeWrite(self.param_ccl) as fd:
                print(f"# Parameter definitions for thorn {self.thornname}",file=fd)
                for name, default, doc, vmin, vmax, options in self.params:
                    t = typeof(default, vmin, vmax)
                    if t == int:
                        print(f'CCTK_INT {name} "{doc}" {{ {vmin}:{vmax} :: "" }} {default}', file=fd)
                    elif t == float:
                        print(f'CCTK_REAL {name} "{doc}" {{ {vmin}:{vmax} :: "" }} {default}', file=fd)
            with SafeWrite(self.interface_ccl) as fd:
                print(f"# Interface definitions for thorn {self.thornname}",file=fd)
                print(f"implements: {self.thornname}",file=fd)
                print(f"inherits: {', '.join(self.external_modules)}",file=fd)
                if gri.ET_driver == "CarpetX":
                    print(f"USES INCLUDE HEADER: loop_device.hxx",file=fd)
                for gf_name in self.gf_names:
                    #TODO: change lines with "if gf_name in ["x","y","z"]:" to check if gftype=="CORE"
                    if gf_name in ["x","y","z"]:
                        continue
                    rhs=gf_name + "_rhs"
                    if re.match(r'.*_rhs',gf_name):
                        tag = "TAGS='checkpoint=\"no\"'"
                    elif rhs in self.gf_names:
                        tag = f"TAGS='rhs=\"{self.thornname}::{rhs}GF\"' "
                    else:
                        tag = ""
                    if self.gf_names[gf_name] == self.thornname:
                        if grid.ET_driver == "Carpet":
                            print(f'REAL {gf_name}GF TYPE=GF TIMELEVELS=3 {{ {gf_name}GF }} "{gf_name}"', file=fd)
                        elif grid.ET_driver == "CarpetX":
                            # CarpetX does not use sub-cycling in time, so it only needs one time level
                            print(f'REAL {gf_name}GF TYPE=GF TIMELEVELS=1 {tag} CENTERING={{ {self.centering[gf_name]} }} {{ {gf_name}GF }} "{gf_name}"', file=fd)

            with SafeWrite(self.schedule_ccl) as fd:
                print(f"# Schedule definitions for thorn {self.thornname}",file=fd)
                for gf_name in self.gf_names:
                    if gf_name in ["x","y","z"]:
                        continue
                    if self.gf_names[gf_name] == self.thornname:
                        if grid.ET_driver == "Carpet":
                            print(f"storage: {gf_name}GF[3]", file=fd)
                        else:
                            print(f"storage: {gf_name}GF[1]", file=fd)
                for src in self.src_files:
                    fsrc = os.path.join(self.src_dir, src.name)

                    readgfs = set()
                    writegfs = set()

                    for func in src.funcs:
                        print(file=fd)
                        # TODO: Add the other special bins here
                        if func.schedule_bin.lower() in ["initial"]:
                            atin = "at"
                        else:
                            atin = "in"
                        print(f"schedule {func.name} {atin} {func.schedule_bin} {{",file=fd)
                        print(f"   LANG: C",file=fd)
                        for readgf in func.readgfs:
                            if readgf in ["x","y","z"]:
                                continue
                            # The symbols in readgfs might not actually be grid
                            # functions. Make sure that they are registered as
                            # such before generating read/write decls.
                            if readgf in self.gf_names:
                                full_name = self.get_full_name(readgf)
                                print(f"   READS: {full_name}(everywhere)",file=fd)
                        for writegf in func.writegfs:
                            if writegf in self.gf_names:
                                full_name = self.get_full_name(writegf)
                                print(f"   WRITES: {full_name}({src.where})",file=fd)
                        if src.where == "interior":
                            for writegf in func.writegfs:
                                if writegf in self.gf_names:
                                    pass #print(f"    SYNC: {full_name}",file=fd)
                        print(f'}} "{func.doc}"',file=fd)
            if not os.path.exists(self.doc_tex):
                doc_src = doc_init
                doc_src = re.sub(r'{thorn}',self.thornname,doc_src)
                doc_src = re.sub(r'{author}',self.author,doc_src)
                doc_src = re.sub(r'{email}',self.email,doc_src)
                doc_src = re.sub(r'{date}',today,doc_src)
                with open(self.doc_tex, "w") as fd:
                    print(doc_src, file=fd)
            with SafeWrite(self.makefile) as fd:
                # for item in outC_function_master_list
                print(makefile_init, end='', file=fd) 
                for src in self.src_files:
                    print(" ",src.name,sep="",end="",file=fd)
                print(file=fd)

            for src in self.src_files:
                fsrc = os.path.join(self.src_dir, src.name)
                # switch to add_to_Cfunction_dict
                #add_to_Cfunction_dict(body='/* body */',includes=['#include <cctk.h>'], \
                #    name='foo',params='/* params */',path_from_rootsrcdir_to_this_Cfunc="/tmp")
                #construct_NRPy_function_prototypes_h("/tmp")
                #construct_NRPy_Cfunctions("/tmp")
                #outCfunction(outfile="x.cc",name='foo2',params='()',body='/* body 2 */')
                with SafeWrite(fsrc) as fd:
                    if gri.ET_driver == "Carpet":
                        print("#include <cctk.h>", file=fd)
                        print("#include <cctk_Arguments.h>", file=fd)
                        print("#include <cctk_Parameters.h>", file=fd)
                    elif gri.ET_driver == "CarpetX":
                        print("#include <fixmath.hxx>", file=fd)
                        print("#include <cctk.h>", file=fd)
                        print("#include <cctk_Arguments.h>", file=fd)
                        print("#include <cctk_Parameters.h>", file=fd)
                        print("#include <loop_device.hxx>", file=fd)
                        print("\nusing namespace Loop;", file=fd)
                    else:
                        assert "Bad value for grid.ET_driver={grid.ET_driver}"
                    for func in src.funcs:
                        print(file=fd)
                        print(f"void {func.name}(CCTK_ARGUMENTS) {{",file=fd)
                        if gri.ET_driver == "Carpet":
                            print(f"  DECLARE_CCTK_ARGUMENTS_{func.name};",file=fd)
                        elif gri.ET_driver == "CarpetX":
                            print(f"  DECLARE_CCTK_ARGUMENTSX_{func.name};",file=fd)
                        print(f"  DECLARE_CCTK_PARAMETERS;",file=fd)
                        print(f"  std::cout << \"Callling {func.name}!!!\" << std::endl;",file=fd)
                        for ii in range(3):
                            print(f"  CCTK_REAL invdx{ii} = 1/CCTK_DELTA_SPACE({ii});",file=fd)
                        print(f"  {func.body}",file=fd)
                        print(f'}}',file=fd)

            # Create the thorn_list entry if needed
            entry = f"{self.arrangement}/{self.thornname}"
            thorn_list = os.path.join("configs", cactus_config, "ThornList")
            self.update_thornlist(entry, thorn_list)
            if cactus_thornlist is not None:
                self.update_thornlist(entry, cactus_thornlist)
        finally:
            if cwd is not None:
                os.chdir(cwd)

    def get_xyz(self):
        assert self.ET_driver == grid.ET_driver
        if grid.ET_driver == "Carpet":
            x,y,z = self.register_gridfunctions("EXTERNAL", ["x","y","z"], external_module="grid")
        elif grid.ET_driver == "CarpetX":
            x,y,z = self.register_gridfunctions("CORE", ["x","y","z"])
        else:
            assert "Bad value for grid.ET_driver={grid.ET_driver}"
        return x,y,z

    def update_thornlist(self, entry, thorn_list):
        if os.path.exists(thorn_list):
            with open(thorn_list, "r") as fd:
                contents = fd.read()
            if re.search(f'^{entry}\\b', contents, re.MULTILINE):
                print(f"Thorn {entry} is already in {os.getcwd()}/{thorn_list}")
            else:
                print(f"Appending {entry} to {os.getcwd()}/{thorn_list}")
                with open(thorn_list, "a") as fd:
                    print(entry, file=fd)
        else:
            print(f"Thornlist {os.getcwd()}/{thorn_list} does not exist. Entry {entry} not added.")
