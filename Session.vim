let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/Documents/Mestrado/DISSERTAÇÃO/tmCOM
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
let s:shortmess_save = &shortmess
if &shortmess =~ 'A'
  set shortmess=aoOA
else
  set shortmess=aoO
endif
badd +77 scripts/sample_runner.py
badd +1 ~/Documents/Mestrado/DISSERTAÇÃO/tmCOM
argglobal
%argdel
$argadd ~/Documents/Mestrado/DISSERTAÇÃO/tmCOM
edit scripts/sample_runner.py
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=99
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
34,40fold
46,48fold
49,51fold
55,57fold
58,60fold
65,68fold
70,71fold
64,72fold
77,78fold
75,79fold
91,95fold
96,103fold
82,110fold
117,119fold
116,119fold
113,127fold
130,131fold
let &fdl = &fdl
let s:l = 77 - ((33 * winheight(0) + 34) / 68)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 77
normal! 0
lcd ~/Documents/Mestrado/DISSERTAÇÃO/thesis
if exists(':tcd') == 2 | tcd ~/Documents/Mestrado/DISSERTAÇÃO/tmCOM | endif
tabnext 1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20
let &shortmess = s:shortmess_save
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
nohlsearch
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
