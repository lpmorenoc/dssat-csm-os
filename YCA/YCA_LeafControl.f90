!!***************************************************************************************************************************
!! This module is intended to calculate behavior of the plant leaf 
!! 10/11/2017 converted from UTF-8 to ANSI
!! Atributes:
!!   
!! Object functions:
!!        
!! Static functions:
!!        
!! Authors
!! @danipilze
!!*********
!
    Module YCA_LeafControl 
    
     USE YCA_First_Trans_m
     USE YCA_Node
     
    contains
    
  
    
    ! true is leaf is active
    logical function isLeafExpanding(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        isLeafExpanding = node%LAGETT <= LLIFGTT
    end function isLeafExpanding
    
    ! true is leaf is active
    logical function isLeafActive(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        isLeafActive = node%LAGETT >= LLIFGTT .AND. node%LAGETT < LLIFGTT+LLIFATT
    end function isLeafActive
    
    ! true is leaf is senescing
    logical function isLeafSenescing(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        isLeafSenescing = node%LAGETT >= LLIFGTT+LLIFATT .AND. node%LAGETT < LLIFGTT+LLIFATT+LLIFSTT
    end function isLeafSenescing
    
    ! true is leaf is alive
    logical function isLeafAlive(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        isLeafAlive = node%LAGETT < LLIFGTT+LLIFATT+LLIFSTT
    end function isLeafAlive
    
    ! true is leaf was active today
    logical function didLeafStartActiveToday(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        didLeafStartActiveToday = node%LAGETT-TTLFLife*EMRGFR  < LLIFGTT .AND. isLeafActive(node)
    end function didLeafStartActiveToday
    
    ! true is leaf started senescing today
    logical function didLeafStartSenescingToday(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        didLeafStartSenescingToday = node%LAGETT-TTLFLife*EMRGFR  < LLIFGTT+LLIFATT .AND. isLeafSenescing(node)
    end function didLeafStartSenescingToday
    
    ! true is leaf will start senescing today
    logical function willLeafStartSenescingToday(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        willLeafStartSenescingToday = isLeafActive(node) .AND. node%LAGETT+TTLFLife*EMRGFR  <= LLIFGTT+LLIFATT                                                      !EQN 371
    end function willLeafStartSenescingToday
    
    ! true is leaf fell today
    logical function didLeafFallToday(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        didLeafFallToday = node%LAGETT-TTLFLife*EMRGFR  < LLIFGTT+LLIFATT+LLIFSTT .AND. .NOT. isLeafAlive(node)
    end function didLeafFallToday
    
    ! real value of the leaf area left to senescence
    real function leafAreaLeftToSenescence(node)
        implicit none
        class (Node_type), intent(in) :: node
        
        leafAreaLeftToSenescence = node%LATL3T - node%LAPS
    end function leafAreaLeftToSenescence
    
    ! set leaf age to active 
    subroutine setLeafAsActive(node)
        implicit none
        class (Node_type), intent(inout) :: node
        
        node%LAGETT = LLIFGTT
    end subroutine setLeafAsActive
    
    ! set leaf age to senescing 
    subroutine setLeafAsSenescing(node)
        implicit none
        class (Node_type), intent(inout) :: node
        
        node%LAGETT = LLIFGTT+LLIFATT                                             !EQN 359
    end subroutine setLeafAsSenescing
    
    ! set leaf age to fall 
    subroutine setLeafAsFall(node)
        implicit none
        class (Node_type), intent(inout) :: node
        
        node%LAGETT = LLIFGTT+LLIFATT+LLIFSTT
    end subroutine setLeafAsFall
    
    ! increase leaf age 
    subroutine leafAge(node)
        implicit none
        class (Node_type), intent(inout) :: node
        
        node%LAGETT = node%LAGETT + TTLFLife*EMRGFR                                              !EQN 358
    end subroutine leafAge

    
    
    END module YCA_LeafControl 