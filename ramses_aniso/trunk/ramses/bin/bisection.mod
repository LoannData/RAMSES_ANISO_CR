  �  G   k820309    4          19.0        �#�]                                                                                                          
       ../amr/bisection.f90 BISECTION                                                     
                                                           
                         @                                '�                   #NGRID    #NPART    #IGRID    #F    #U    #FP 	   #UP 
                �                                                               �                                                             �                                                                       &                                                       �                                          P                             &                   &                                                       �                                         �                 
            &                   &                                                       �                             	                                        &                   &                                                       �                             
            p                
            &                   &                                                                                                                                                            @                                                      @                                                                    &                   &                                                    @                                                   
                &                                                                                                                                              3         @                                                                    &                                                                                                                                                                                                         
                                                                                                                                                                                                                                                 @                                                                    &                   &                                                                                                                                                                                                                                                                      
                                                       
                @                                                   
                &                   &                                                    @                                                   
                &                   &                                                    @                                                    
                &                   &                                                    @                                !                   
                &                   &                                                    @ @                               "                                   &                                                                                       #                                                        $                     @                                 %                                   &                                                                                       &                                                        '                     @                                 (                                   &                                                    @                                 )                                   &                                                    @                                 *                                   &                                                                                        +                                                        ,                                                      -                   
                &                   &                                                                                        .                                                       32                                             /                                                                                                     0                                                       1                                   &                                                                                      2                                   &                                                    @                                 3                                   &                                                                                       4                                                       5            �                       &                                           #COMMUNICATOR                                               6            �                       &                   &                                           #COMMUNICATOR    #         @                                   7                    #X 8   #C 9   #NN :             
                                 8                   
              &                   &                                                     D                                 9                                  &                                                     
                                  :           #         @                                   ;                    #UPDATE <             
                                  <           #         @                                  =                     #         @                                  >                    #LEV ?   #DIR @             
                                  ?                     
                                  @           %         @                                A                    
       #X B             
                                 B     
      #         @                                  C                    #LEV D   #DIR E   #WALLS F             
                                  D                     
                                  E                     
                                 F                   
              &                                              �   '      fn#fn    �   @   j   AMR_PARAMETERS      @   J   AMR_COMMONS )   G  �       COMMUNICATOR+AMR_COMMONS /   �  H   a   COMMUNICATOR%NGRID+AMR_COMMONS /     H   a   COMMUNICATOR%NPART+AMR_COMMONS /   f  �   a   COMMUNICATOR%IGRID+AMR_COMMONS +   �  �   a   COMMUNICATOR%F+AMR_COMMONS +   �  �   a   COMMUNICATOR%U+AMR_COMMONS ,   R  �   a   COMMUNICATOR%FP+AMR_COMMONS ,   �  �   a   COMMUNICATOR%UP+AMR_COMMONS "   �  p       DP+AMR_PARAMETERS '     @       BISEC_ROOT+AMR_COMMONS '   Z  �       BISEC_NEXT+AMR_COMMONS '   �  �       BISEC_WALL+AMR_COMMONS $   �  q       NDIM+AMR_PARAMETERS '   �  �       BISEC_INDX+AMR_COMMONS (   �  @       NBINODES+AMR_PARAMETERS ,   �  @       NBILEAFNODES+AMR_PARAMETERS &   	  @       BOXLEN+AMR_PARAMETERS +   G	  @       ICOARSE_MAX+AMR_PARAMETERS +   �	  @       ICOARSE_MIN+AMR_PARAMETERS '   �	  @       VERBOSE+AMR_PARAMETERS !   
  @       NCPU+AMR_COMMONS '   G
  �       BISEC_HIST+AMR_COMMONS '   �
  @       BISEC_NRES+AMR_COMMONS +   +  @       NBILEVELMAX+AMR_PARAMETERS !   k  @       MYID+AMR_COMMONS &   �  @       BISEC_RES+AMR_COMMONS )   �  @       BISEC_TOL+AMR_PARAMETERS .   +  �       BISEC_CPUBOX_MIN2+AMR_COMMONS .   �  �       BISEC_CPUBOX_MAX2+AMR_COMMONS -   s  �       BISEC_CPUBOX_MIN+AMR_COMMONS -     �       BISEC_CPUBOX_MAX+AMR_COMMONS +   �  �       BISEC_CPU_LOAD+AMR_COMMONS +   G  @       JCOARSE_MIN+AMR_PARAMETERS +   �  @       KCOARSE_MIN+AMR_PARAMETERS ,   �  �       NEW_HIST_BOUNDS+AMR_COMMONS "   S  @       NX+AMR_PARAMETERS "   �  @       NY+AMR_PARAMETERS .   �  �       BISEC_HIST_BOUNDS+AMR_COMMONS +   _  �       BISEC_IND_CELL+AMR_COMMONS '   �  �       CELL_LEVEL+AMR_COMMONS $   w  @       NCOARSE+AMR_COMMONS (   �  @       NGRIDMAX+AMR_PARAMETERS    �  �       XG+AMR_COMMONS '   �  r       NVECTOR+AMR_PARAMETERS )     p       TWOTONDIM+AMR_PARAMETERS "   }  @       NZ+AMR_PARAMETERS $   �  �       CPU_MAP+AMR_COMMONS     I  �       SON+AMR_COMMONS "   �  �       FLAG1+AMR_COMMONS )   a  @       NLEVELMAX+AMR_PARAMETERS #   �  �       ACTIVE+AMR_COMMONS &   ?  �       RECEPTION+AMR_COMMONS %   �  ^       CMP_BISECTION_CPUMAP '   S  �   a   CMP_BISECTION_CPUMAP%X '   �  �   a   CMP_BISECTION_CPUMAP%C (   �  @   a   CMP_BISECTION_CPUMAP%NN     �  T       BUILD_BISECTION '     @   a   BUILD_BISECTION%UPDATE )   W  H       INIT_BISECTION_HISTOGRAM *   �  Z       BUILD_BISECTION_HISTOGRAM .   �  @   a   BUILD_BISECTION_HISTOGRAM%LEV .   9  @   a   BUILD_BISECTION_HISTOGRAM%DIR #   y  W       ROUND_TO_BISEC_RES %   �  @   a   ROUND_TO_BISEC_RES%X .     e       SPLITSORT_BISECTION_HISTOGRAM 2   u  @   a   SPLITSORT_BISECTION_HISTOGRAM%LEV 2   �  @   a   SPLITSORT_BISECTION_HISTOGRAM%DIR 4   �  �   a   SPLITSORT_BISECTION_HISTOGRAM%WALLS 