function [S,S_p,D_xp,D_xn,pools,S_x] = Stoi(r)

    w = r / (1 - r);
        %1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	
   S = [-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	w	0	0	0	0	0	0	0	0	0	0	0	0	0  %1
        1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	-w	0	0	0	0	0	0	0	0	0	0	0	0  %2
        0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	w	0	0	0	0	0	0	0	0	0	0	0  %3        
        0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0	0	0	0	0	0	0	0	0	0  %4        
        0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	w	0	0	0	0	0	0	0	0	0  %5        
        0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	w	0	0	0	0	0	0	0	0  %6       
        0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0	0	0	0	0	0	0  %7        
        0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0	0	0	0	0	0  %8        
        0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  %9        
        0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0  %10        
        0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  %11        
        0	0	0	0	0	1	-1	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0  %12        
        0	0	0	0	0	0	1	-1	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0  %13        
        0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0  %14        
        0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  %15        
        0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0  %16        
        0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0  %17        
        0	0	0	0	0	0	0	0	0	0	1	0	0	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0  %18        
        0	0	0	0	0	0	0	0	0	0	0	0	1	1	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0  %19        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0  %20        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1  %21        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0  %22        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0  %23        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0  %24        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0  %25        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  %26        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  %27        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0	0	0	0	0  %28        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0	0	0	0  %29	       
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  %30
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0	0	0  %31        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	-1	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	-w	0	0  %32
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w	0  %33        
        0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-w];%34

          %1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32  33  34
    S_p = [1	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	1	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	1	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	1	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	1	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	1	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	1	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	1   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0           
           0	0	0	0	0	0	0	0	0	0	0	1	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	1	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	1	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	1	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	1	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0   0   0   0   
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	1	0	0	0	0	0	0	0   0   0   0   0           
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	1	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	1	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	1	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	1	0	0	0   0   0   0   0
           0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	1	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	1   0   0   0   0   
           0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   1   0   0   0           
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0   0   0   0   
           0	0	0	1	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   1   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   1   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   1   0   0           
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   1
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   1
           0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0           
           0	0	0	0	0	0	0	0	0	0	0	1	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	1	0	0	0	0	0	0	0	0	0	0   0   0   0   0           
           0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
           0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   1   0   0];

           %1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32  33  34
    D_xp = [0	0	0	0	0	0	0	0	0	1	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   1	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	1	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0            
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	1	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	1	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	1	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	1	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0            
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	1	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   1   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   1   0   0            
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   1   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   1];
        
           %1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32  33  34
    D_xn = [1	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	1	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	1	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0            
            0	0	0	0	1	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	1	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	1	0	0	0	0	0   0   0   0   0            
            0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	1	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	1	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	1	0	0	0	0	0	0	0	0	0	0   0   0   0   0            
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	1	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
            0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0   0   0   0];

            %1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32  33  34
    pools = [0	0	0	0	0	0	0	0	r	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
             1-r 0	0	0	0	0	0	0	0	r	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
             0	0	0	0	0	0	0	0	0	0	r	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
             0	1-r	0	0	0	0	0	0	0	0	0	r	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0             
             0	0	0	0	0	0	0	0	0	0	0	0	r	0	0   0	0	0	0	0	0	0	0	0	0	0	0	1-r	0	0   0   0   0   0             
             0	0	0	0	0	0	0	0	0	0	0	0	0	r	0   0	0	0	0	0	0	0	0	0	0	0	0	0	1-r	0   0   0   0   0
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	r   0	0	0	0	0	0	0	0	0	0	0	0	0	0	1-r 0   0   0   0             
             0	0	1-r	0	0	0	0	0	0	0	0	0	0	0	0   r	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0
             0	0	0	1-r	0	0	0	0	0	0	0	0	0	0	0   0	r	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0             
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	r	0	0	0	0	0	0	0	0	0	0	0	0   1-r 0   0   0
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	r	0	0	0	0	0	0	0	0	0	0	0   0   0   1-r 0             
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	r	0	0	0	0	0	0	0	0	0	0   0   1-r 0   0
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	r	0	0	0	0	0	0	0	0	0   0   0   0   1-r             
             0	0	0	0	1-r	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	r	0	0	0	0	0	0	0	0   0   0   0   0
             0	0	0	0	0	1-r	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	r	0	0	0	0	0	0	0   0   0   0   0             
             0	0	0	0	0	0	1-r	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	r	0	0	0	0	0	0   0   0   0   0
             0	0	0	0	0	0	0	1-r	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	r	0	0	0	0	0   0   0   0   0             
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	r	0	0	0	0   0   0   0   0
             0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	r	0	0	0   0   0   0   0];
    
          %1   2   3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56
    S_x = [1   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	1	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   1   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0 	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   1   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   1	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	1	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0           
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	1	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	1	0	0	0	0	0   0           
           0   0   0	0	1	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	1	0	0	0	0   0	0	0	0	0	0	0	0   0           
           0   0   0	0	0	1	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	1	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	1	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	1	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	1	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0           
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	1	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	1	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0           
           0   0   0	0	0	0	0	0	0	0	0	0	0	1	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   1	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	1	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	1	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0         
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	1   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	1	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	1	0	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	1	0	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	1	0	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	1	0	0	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0          
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	1	0	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	1	0	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	1	0	0   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0           
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	1   0   0   0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   1	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   1   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0
           0   0   0	0	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0	0	0	0	0	0	0	0   0   0   0   0   0	1	0	0	0	0	0	0	0	0	0	0	0   0	0	0	0	0	0	0	0   0];  
      
end
