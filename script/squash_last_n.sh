#!/bin/bash

# 合并最近 N 个提交
# 用法： ./squash_last_n.sh N

# if [ $# -ne 1 ]; then
#     echo "用法: $0 <number-of-commits>"
#     exit 1
# fi

# N=$1

# # 创建一个临时文件
# TMPFILE=$(mktemp)

# # 生成 rebase 脚本
# echo "生成 rebase 指令..."

# git log -n $N --pretty=format:"pick %H %s" | tac > $TMPFILE

# # 除第一行外，把所有 pick 改成 squash
# sed -i '2,$ s/^pick /squash /' $TMPFILE

# echo "下面是生成的 rebase 脚本："
# cat $TMPFILE

# echo
# echo "准备开始 rebase..."

# # 启动交互式 rebase
# GIT_SEQUENCE_EDITOR="cat $TMPFILE >" git rebase -i HEAD~$N

# echo "完成后如有需要，可执行：git push origin 分支名 --force"
git reset --soft HEAD~1
git commit --amend
git push --force