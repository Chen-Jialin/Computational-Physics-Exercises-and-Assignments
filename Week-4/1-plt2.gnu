plot "1-pdf.txt" u 1:i w l lw 3 title sprintf("t = %f",i * 0.02)
i=i+1 #变量加一
pause 0.05
if(i<100) reread