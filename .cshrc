set path=($path . ~/Documents/bin /home/gworthey/Documents/bin \
  /home/gworthey/Documents/data/misc/iraf2/imfort)

setenv OMP_NUM_THREADS 4
setenv MESASDK_ROOT ~/Documents/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.csh

if ( $?prompt ) then
	set history=32
endif

setenv iraf /home/gworthey/iraf/
source /home/gworthey/iraf/unix/hlib/irafuser.csh

#/usr/c/kurucz/synthe/atlas/bin /usr/c/kurucz/synthe/atlas/src/CDROM/bin \

alias work 'cd /home/student/Documents/mesa3723/star'
limit coredumpsize 0

alias ls 'ls -aF'
alias lst 'ls -alF | more'
alias dls 'ls -aF | grep /'
alias rm 'rm -i'
alias cd 'cd \!*; prompt'
alias prompt  'set prompt="`hostname` $cwd>"'
prompt
alias wipe 'rm .*~ *~ *.log *.bak *.trace fort.* *.dvi *.aux'
alias dir 'ls -aF'
alias smemacs 'emacs -fn 6x12'
alias lemacs  'emacs -fn 9x15bold'
alias ffind   'find . -name "\!*" -print'
#  cat sample.pdf | acroread -toPostScript > sample.ps 
#  or acroread -toPostScript -level2 pdf_file_1
#alias pdfps   'acroread -toPostScript -level2 '
#alias xterm 'xterm -bg white'
#alias xgterm 'xgterm -bg white'

alias mesa 'cd ~/Documents/mesa7624/mesa/star'
alias tong 'cd ~/Documents/tang'

#setenv KURUCZ_LIB /usr/c/kurucz/synthe/atlas/lib

#set LD_LIBRARY_PATH=( /usr/local/lib /lib /usr/lib ./ /usr/X11R6/lib )

setenv PGPLOT_DIR /usr/lib/pgplot5
setenv PGPLOT_FONT /usr/lib/pgplot5/grfont.dat
setenv SPS_HOME /home/student/fsps

setenv oref /home/student/Documents/hst/raw/

alias ur_setup 'eval `/home/gworthey/.ureka/ur_setup -csh \!*`'
alias ur_forget 'eval `/home/gworthey/.ureka/ur_forget -csh \!*`'

ur_setup

