%%%%%%%%%%%%%%%
%%This is Jeffrey's personal article class
%%For personal use only
%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{JeffArticle}[2019/03/03]
%This is optional. It's only used to check the version of class

\LoadClass{article}
\RequirePackage{titling}
\RequirePackage{xpatch}
\RequirePackage{graphicx}

%%控制标题格式
\RequirePackage{titlesec}
\titleformat{\section}[hang]
{\Large\scshape\raggedright}
{\itshape \thesection\quad}
{0em}{}[\titlerule]
\titleformat{\subsection}{\large \raggedright}{$\bullet$ {\large \itshape \thesubsection \ } }{0em}{}

%\makeatletter
%%控制页面格式
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\lhead{\@title}
\chead{}
\rhead{\thepage}
\lfoot{jeffreyyao@pku.edu.cn}
\cfoot{}
\rfoot{\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0.4pt}



\renewcommand\maketitle{%
	%\hfill
	
	\begin{minipage}{0.95\textwidth}
		\vskip 2em
		%\rule{\linewidth}{0.4pt}
		{\LARGE \centering \textbf \@title \par }
		%\rule{\linewidth}{0.4pt}
		\vskip 1 em
		{\large \centering \@author \par}
		{\large \centering \@date \par}
	\end{minipage}
	\vskip 1em \par
}
\newcommand{\titlepagestyle}{fancy}
\xpatchcmd{\maketitle}{\thispagestyle{plain}}{\thispagestyle{\titlepagestyle}}{}{}  %让首页加上fancy效果  默认为plain
