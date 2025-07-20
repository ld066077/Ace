
- **MPI 进程数**由 `-n $Ntasks` 决定。
- **每个 MPI 进程的 OpenMP 线程数**由 `-c $CoresPerTask`（或者 `OMP_NUM_THREADS` 环境变量）决定。
```bash
srun -N $Nnodes -n $Ntasks -c $CoresPerTask -p debug3N ./NLT >& ./run.out
```
##### 正确使用例子
主程序
```fortran
program main
    use global_data
    use mpi
    
    use numerical_module
    implicit none
    integer :: times, k, j, col
    integer :: rank, size, ierr
    real(8) :: min_errN, flag_omega, omega
    real(8) :: send_buffer, recv_buffer
    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    call initialize_globals()
    ! Set initial conditions for omega
    omega = 0.0
    times = int((omega_end - omega_start) / (0.01))
    ! Main loop to perform shooting method for omega
    if (rank /= 0) then
        do k = rank, times, size-1
            omega = omega_end - (k - 1) * domega
            flag_omega = omega
            call find_minimum_error(flag_omega)
            send_buffer = flag_omega
            call MPI_Send(send_buffer, 1, MPI_DOUBLE_PRECISION, 0, k, MPI_COMM_WORLD, ierr)
        end do
    else
        ! Rank 0 receives and processes results
        open(unit=10, file='mineData.dat', status='unknown', position='append')
        do k = 1, times
            call MPI_Recv(recv_buffer, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            ! Convert integer 'k' to a real(8) before writing
            write(10, '(F25.10, 1X, ES25.10E4)') recv_buffer, real(k, 8)
                print *, 'Received valid iteration k = ', k, '/', times
        end do
        close(10)
    endif
    call MPI_Finalize(ierr)

end program main
```
**MPI_TAG**用k，表示进程序号，输出结果按照k值顺序输出，发送接收无误，顺序正确
调用模块
```Fortran
module global_data
    implicit none
    real(8) :: domega, omega_start, omega_end
contains
    ! Initialization subroutine to set up all global variables
    subroutine initialize_globals()
       ! integer :: i
        domega = 0.01d0
        omega_start = 1.0
        omega_end = 2.0
    end subroutine initialize_globals
end module global_data
```

```
module numerical_module
    use global_data
    implicit none
    contains
    subroutine find_minimum_error(flag_omega)
        implicit none
        real(8), intent(inout) :: flag_omega
        flag_omega = flag_omega + domega
    end subroutine find_minimum_error
  end module numerical_module
```
编译运行命令：
```
gfortran -g -c MPI_YZGLOBAL.f90
gfortran -g -c MPI_YZNUM.f90
mpif90 -g -c MPI_YZ.f90
mpif90 -g -o MPI_YZ MPI_YZGLOBAL.o MPI_YZNUM.o MPI_YZ.o /lib64/liblapack.so.3 /lib64/libblas.so.3

sun ./MPI_YZ
```
##### 基本概念
**并行计算机分类**
**SIMD** 单指令多数据并行计算机  如对单个数组操作，其所有元素同时处理
**MIMD** 多指令多数据
**SPMD** 单程序多数据
**MPMD** 多程序多数据

**并行编程模型**
数据并行
消息传递
...

	对于机群计算，单次通信时间>>单次计算，降低通信次数，以计算换通信。理想情况：计算与通信的重叠
MPI 并行函数库，使得串行程序如Fortran、C调用后扩展为并行程序。为**消息传递编程模型**，服务于进程通信。 
##### MPI通信模式
```
在 MPI 程序中，特别是在点对点通信（如 `MPI_Send` 和 `MPI_Recv`）中，存在一种机制可以避免因进程间通信的不同步而引起的问题。这种机制涉及到 MPI 的通信模式（阻塞和非阻塞通信）的使用。

### 问题的理解

你提到的疑问是关于 MPI 程序中的一个常见问题：如果接收进程（进程 0）调用 `MPI_Recv` 早于其他进程（非 0 进程）调用 `MPI_Send`，会发生什么情况？

在这个特定的程序中，进程 0 使用 `MPI_Recv` 从其他进程接收数据，而其他进程使用 `MPI_Send` 发送数据到进程 0。那么，问题就是如果接收操作先执行，发送操作还没有准备好，程序会如何处理？

### MPI 的通信机制

在 MPI 的标准通信模式中，`MPI_Send` 和 `MPI_Recv` 都是阻塞的操作。这意味着：

1. **`MPI_Recv`**（阻塞接收）：
   - 当进程 0 调用 `MPI_Recv` 时，如果消息还没有准备好，它会阻塞（暂停）直到消息到达。这意味着进程 0 会等待，直到有一个符合接收条件的消息被发送过来。
   
2. **`MPI_Send`**（阻塞发送）：
   - 当非 0 进程调用 `MPI_Send` 时，它会阻塞直到消息可以被发送。通常，这意味着它会等待直到 MPI 系统确认接收进程有足够的空间来接收这条消息，或者直到消息被缓冲（如果有足够的缓冲空间）。

### 为什么不会出问题

在 MPI 中，发送和接收的阻塞行为保证了即使接收命令先于发送命令，程序也不会出错。具体来说：

- 如果进程 0 调用 `MPI_Recv` 早于其他进程的 `MPI_Send`，进程 0 会阻塞，等待消息到达。
- 当其他进程后续调用 `MPI_Send` 时，这些消息将会发送到进程 0，并解除 `MPI_Recv` 的阻塞状态，继续执行。

因此，在这种阻塞模型下，无论进程 0 的 `MPI_Recv` 调用是否早于其他进程的 `MPI_Send`，程序都会正常运行。

### 非阻塞通信的考虑

虽然这个程序在阻塞通信模式下运行正常，但在一些情况下，程序可能需要提高效率或者避免长时间阻塞。在这种情况下，MPI 提供了非阻塞通信（如 `MPI_Isend` 和 `MPI_Irecv`），这允许程序在消息发送或接收操作进行的同时继续执行其他操作。然而，这需要更加仔细的管理，确保发送和接收的正确匹配，以及在适当的时候完成这些非阻塞操作（通过 `MPI_Wait` 或 `MPI_Test` 等函数）。

### 总结

因此，尽管接收进程 0 的 `MPI_Recv` 可能先于其他进程的 `MPI_Send`，在标准阻塞模式下，MPI 的行为保证了程序的正确性和同步性。程序会等待直到数据可用，从而避免了通信上的竞争条件和数据丢失的问题。
```



高精度的计算可能在某些表达式中产生过大的中间值（溢出）或过小的中间值（下溢），影响数值稳定性。
例如较高精度时LU分解出错，低精度没有出错，精度过高确实可能导致中间计算结果出现极小值甚至接近零的情况，从而引发数值不稳定问题，尤其是在对矩阵进行分解（如 LU 分解）时。这种情况通常由浮点数下溢或过小的计算值对矩阵条件数的影响引起。
**通过减少 `s` 和 `alpha` 的精度（如保留前四位有效数字），可以适当增加中间值的误差范围，使计算结果不至于过于接近零。**
- 矩阵的条件数（condition number）是数值稳定性的关键指标。当矩阵的元素分布范围过于极端（如包含极小或极大的数值），条件数会变大，数值稳定性下降。
- 降低 `s` 和 `alpha` 的精度会减小计算中的数值差异，使矩阵的条件数减小，从而提高 LU 分解的成功率。
- 如果问题的核心需求是**高精度结果（如科学研究）**，建议在**条件数优化和矩阵扰动**的基础上再**提高数值稳定性**，**而非简单降低精度**。
- 1. **动态选择适当的精度**：
    
    - 如果可以接受较大的误差，减少精度是**简单**有效的方法。
    - 对于非常敏感的系统，需**平衡计算稳定性和结果精度**。
在 LU 分解中添加条件数检查：
```fortran
! Compute condition number and check
condition_number = calculate_condition_number(CoeMatrix)
if (condition_number > 1.0e12_dp) then
    print *, "Warning: Condition number is too large:", condition_number
    stop
end if

function calculate_condition_number(A) result(cond)
    implicit none
    complex(dp), intent(in) :: A(:,:)
    real(dp) :: cond
    real(dp), allocatable :: singular_values(:)
    call zgesvd('N', 'N', size(A, 1), size(A, 2), A, size(A, 1), singular_values, ...)
    cond = maxval(singular_values) / minval(singular_values)
end function calculate_condition_number


```



```bash
曙光集群:
gfortran -g -c global_data.f90
gfortran -g -c FTLR_module.f90
gfortran -g -c numerical_module.f90
gfortran -g -c shoot_omega.f90                                                                \                                        mpif90 -g -c shoot_omega.f90
gfortran -g -o shoot_omega global_data.o FTLR_module.o numerical_module.o shoot_omega.o /lib64/liblapack.so.3 /lib64/libblas.so.3

sbatch submit_shoot_nompi_omega.sh

gdb ./shoot_omega
break numerical_module.f90:10

run




gfortran -fopenmp -g -o shoot_omega global_data.o FTLR_module.o numerical_module.o shoot_omega.o /lib64/liblapack.so.3 /lib64/libblas.so.3
2
srun


mpif90 -g -o shoot_omega global_data.o FTLR_module.o numerical_module.o shoot_omega.o /lib64/liblapack.so.3 /lib64/libblas.so.3

ubuntu:
gfortran -g -c global_data.f90
gfortran -g -c FTLR_module.f90
gfortran -g -c numerical_module.f90
gfortran -g -c phi_2.f90
gfortran -g -o phi_2 global_data.o FTLR_module.o numerical_module.o phi_2.o /lib/x86_64-linux-gnu/liblapack.so.3 /lib/x86_64-linux-gnu/libblas.so.3
```


#### bash&sbatch 脚本

[script_bash_sbatch](file:///F:%5Ccode%5Cscript_bash_sbatch)
##### sugon
sub.sh
```bash
Name=GAM_KE
Partition=hfacnormal01
Nnodes=2
CoresPerNode=32
#Queue=TH_HPC3N
#Queue=debug3N
if [ -f "input/mesh.dat" ]; then
  Ntasks=$(awk '/nmu/{gsub(/[^0-9]/,"",$0);print $0}' input/mesh.dat)
  Nionnum=$(awk '/IonNum/{gsub(/[^0-9]/,"",$0);print $0}' input/mesh.dat)
else
  echo "input/mesh.dat is not exist!"
  exit 1
fi
TasksPerNode=$(($Ntasks/$Nnodes))
CoresPerTask=$(($CoresPerNode/$TasksPerNode))
sed -i "s/#SBATCH --job-name=.*/#SBATCH --job-name=$Name/" job.sh
sed -i "s/#SBATCH --partition=.*/#SBATCH --partition=$Partition/" job.sh
sed -i "s/#SBATCH --nodes=.*/#SBATCH --nodes=$Nnodes/" job.sh
sed -i "s/#SBATCH --ntasks-per-node=.*/#SBATCH --ntasks-per-node=$TasksPerNode/" job.sh
sed -i "s/#SBATCH --cpus-per-task=.*/#SBATCH --cpus-per-task=$CoresPerTask/" job.sh
sed -i "s/Ntasks=.*/Ntasks=$Ntasks/" job.sh
sed -i "s/Nionnum=.*/Nionnum=$Nionnum/" job.sh
sed -i "s/Nnodes=.*/Nnodes=$Nnodes/" job.sh
sed -i "s/CoresPerNode=.*/CoresPerNode=$CoresPerNode/" job.sh
sed -i "s/TasksPerNode=.*/TasksPerNode=$TasksPerNode/" job.sh
sed -i "s/CoresPerTask=.*/CoresPerTask=$CoresPerTask/" job.sh
HOME=$(pwd)
if [ ! -d "${HOME}/dump" ]; then
  mkdir ${HOME}/dump
fi
if [ ! -d "${HOME}/F0file" ]; then
  mkdir ${HOME}/F0file
fi
sbatch job.sh
```
job.sh
```sbatch
#!/bin/sh
#SBATCH --job-name=GAM_KE
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --partition=hfacnormal01
DIREC=${HOME}/${SLURM_JOB_ID: -5}_$SLURM_JOB_NAME
mkdir $DIREC
ulimit -s unlimited
cd $DIREC
cp -r ${HOME}/NLT ./
cp -r ${HOME}/input ./
mkdir srcs
cp ${HOME}/*.f90 srcs/
Ntasks=32
Nionnum=1
Nnodes=2
CoresPerNode=32
TasksPerNode=16
CoresPerTask=2
echo "Total Tasks: $Ntasks"
echo "Total Nodes: $Nnodes"
echo "Cores per Node: $CoresPerNode"
echo "Tasks per Node: $TasksPerNode"
echo "Cores per Task: $CoresPerTask"
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
echo The job starts at time :
date
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
module purge
module load mpi/intelmpi/2017.4.239
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export OMP_NUM_THREADS=2 #设置openMP使用的线程数
srun ./NLT >& ./run.out
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
echo The job ends at time :
date
echo "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-"
```


## 为什么单核串行运行就不报错？ 后续操作：仅0rank计算，不需要Bcast
报错：
```text
程序调试： 1 # Intel ifort / mpiifort: 
mpiifort -check bounds -traceback -g -O0 main.f90

2 如果您的编译器支持浮动异常处理，请打开它（例如，Intel Fortran 中的`-fpe0`，当您执行“除以零”或产生 Inf 时，这通常会立即停止

3 您有通过调用`WOFZ`返回的`flag` 。如果`flag`为`.true.`如果超出范围，至少打印警告或进行检查。现在，您永远不会检查`flag`

segmentation faults (SIGSEGV)段错误 可能原因： 数组越界，除0，矩阵奇异，
```
很多时候，一个**本来就存在**但被隐藏的越界/非法访问，在单线程下可能“碰巧”访问到尚未被其他数据占用的内存区域，程序也就没有立刻崩溃。并行后：

1. **内存布局**发生了变化  
    并行编译或进程/线程启动后，分配内存时对变量的排布顺序和地段往往和串行不一样；越界写就更可能破坏到**关键区域**(比如并行库的内部数据结构)，因此更容易触发段错误。
    
2. **时间顺序**变化  
    多线程下可能有多条线程/进程几乎同时在更新/读取某段内存。若代码里存在对同一个数组的错误操作或缺失同步，就可能导致段错误。
    
3. **编译优化**差异  
    并行编译常常带有更激进的优化(如`-O2/-O3`)，而单核调试版本可能是`-O0`或者带 bounds-check 之类保护；在高优化下，某些“未定义行为”更易暴露出崩溃。
    

因此，“单核正常、并行崩溃”并不一定说明是因为“并行访问量过大导致内存不够”，而是往往暗示你**原本就存在的漏洞**(比如越界访问、未初始化变量等)在并行环境下被放大或更快触发。




| ==111== | ==121== | ==131== |
| :---------: | :---------: | :---------: |
|      211      |     221      |     231      |
|     311      |     321      |     331      |
|     ==112==      |     ==122==      |     ==132==      |
|     212      |     222      |     232      |
|     312      |     322      |     332      |
|     ==113==      |     ==123==      |     ==133==      |
|     213      |     223      |     233      |
| 313 | 323 | 333 |
Fortran 矩阵存储顺序：第一层第一列，第二列，第三列。第二层第一列......
MPI分解rank 放在最后一维