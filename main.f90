PROGRAM MAIN
COMPLEX :: roots(3)
EXTERNAL fn
roots = 0
CALL MULLER(fn,roots,3,200,0.1e-5,0.1e-6)
PRINT*, (roots(i), i=1,3)
END PROGRAM MAIN
