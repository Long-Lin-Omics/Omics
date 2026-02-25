#!/bin/bash

# 用法检查
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <script_file> <lines_per_task> <cpus_per_task>"
    exit 1
fi

SCRIPT_FILE="$1"
LINES_PER_TASK="$2"
CPUS_PER_TASK="$3"

# 前缀
SLURM_PREFIX="slurm_job"

# 检查脚本文件是否存在
if [ ! -f "$SCRIPT_FILE" ]; then
    echo "Error: file '$SCRIPT_FILE' not found"
    exit 1
fi

# 初始化
task_num=1
cmd_buffer=""
line_count=0

while IFS= read -r line; do
    # 拼接命令
    if [ -z "$cmd_buffer" ]; then
        cmd_buffer="$line"
    else
        cmd_buffer="$cmd_buffer"$'\n'"$line"
    fi
    ((line_count++))

    # 到达每批行数，提交 Slurm
    if [ "$line_count" -eq "$LINES_PER_TASK" ]; then
        SLURM_FILE="${SLURM_PREFIX}_${task_num}.slurm"
        cat > "$SLURM_FILE" <<EOL
#!/bin/bash
#SBATCH --job-name=task_$task_num
#SBATCH --output=task_$task_num.out
#SBATCH --error=task_$task_num.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --partition=highmem

$cmd_buffer
EOL

        sbatch "$SLURM_FILE"

        # 重置
        cmd_buffer=""
        line_count=0
        ((task_num++))
    fi
done < "$SCRIPT_FILE"

# 提交最后一批（如果命令数不能整除每批行数）
if [ -n "$cmd_buffer" ]; then
    SLURM_FILE="${SLURM_PREFIX}_${task_num}.slurm"
    cat > "$SLURM_FILE" <<EOL
#!/bin/bash
#SBATCH --job-name=task_$task_num
#SBATCH --output=task_$task_num.out
#SBATCH --error=task_$task_num.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$CPUS_PER_TASK
#SBATCH --time=01:00:00
#SBATCH --partition=standard

$cmd_buffer
EOL

    sbatch "$SLURM_FILE"
fi