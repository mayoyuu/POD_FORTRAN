#!/bin/bash

# ====================================================
# 项目重构脚本: 将 CAT_Fortran 架构重命名为 POD
# 功能: 批量替换文件内容、重命名文件、清理编译缓存
# ====================================================

echo "🚀 开始将项目前缀从 'cat_' 重构为 'pod_'..."

# 1. 替换文件内容中的文本
echo "📝 正在扫描并替换 .f90 和 fpm.toml 文件内容..."
# 考虑到可能存在大写的情况，同时替换小写和大写
# 这里的代码兼容了 Linux (GNU sed) 和 macOS (BSD sed)
if sed --version >/dev/null 2>&1; then
    # Linux / WSL 环境
    find . -type f \( -name "*.f90" -o -name "fpm.toml" \) -exec sed -i 's/cat_/pod_/g' {} +
    find . -type f \( -name "*.f90" -o -name "fpm.toml" \) -exec sed -i 's/CAT_/POD_/g' {} +
else
    # macOS 环境
    find . -type f \( -name "*.f90" -o -name "fpm.toml" \) -exec sed -i '' 's/cat_/pod_/g' {} +
    find . -type f \( -name "*.f90" -o -name "fpm.toml" \) -exec sed -i '' 's/CAT_/POD_/g' {} +
fi
echo "✅ 文件内容替换完成！"

# 2. 重命名包含 'cat_' 的文件或文件夹
echo "📁 正在重命名物理文件..."
# 使用 -depth 确保自下而上处理（先改子文件，再改父文件夹，防止路径失效）
find . -depth -name '*cat_*' | while read -r filepath; do
    dir=$(dirname "$filepath")
    base=$(basename "$filepath")
    # 字符串替换: cat_ -> pod_
    newbase="${base//cat_/pod_}"
    mv "$filepath" "$dir/$newbase"
    echo "   [已重命名] $base  ->  $newbase"
done
echo "✅ 文件重命名完成！"

# 3. 清理旧的编译缓存
echo "🧹 正在清理旧的 fpm build 缓存..."
if [ -d "build" ]; then
    rm -rf build/
    echo "✅ 缓存已清理！"
fi

echo "===================================================="
echo "🎉 重构大功告成！现在你的项目是 POD 了。"
echo "💡 建议运行 'fpm build' 来验证编译是否顺利。"