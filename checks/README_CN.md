完成编译后，从包根目录执行：

```sh
cp -R checks/static_equilibrium run_static
NATIVE_WIRING_CHECKS=1 MPI_OVERSUBSCRIBE=1 bash run.sh run_static
```

这是原配置静止显式一步（8rank），比例校正关闭。正式计算不需要NATIVE_WIRING_CHECKS。完整研究分析脚本与原始数组留在原研究目录。
