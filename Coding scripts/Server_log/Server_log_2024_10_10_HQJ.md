### **Date**
>20241010
### **Environment**
>server:10.30.62.1
### **Task**
* Mount net drive: 浪潮服务器10.120.17.34
* 修改AVITI用户名密码

### **Log**
```shell
# 尝试mount浪潮服务器10.120.17.34
sudo mount -t cifs //10.120.17.34/BioCRF /admin_file -o username=BioCRF,password=biocrf@202401,vers=3.0
# 成功

# 尝试设置开机自动挂载
vim /etc/fstab
//10.120.17.34/BioCRF /admin_file cifs username=BioCRF,password=biocrf@202401,vers=3.0
# Ctrl + X, Y, Enter to save edit
```

### **Environment**
>server:AVITI
### **Task**
* 重命名用户guochao为biocrfgz并修改密码为biocrfgz

### **Log**
```shell
su - # switch to root
kill -u guochao
sudo usermod -l biocrfgz guochao
# CLI界面 → 设置 → 用户管理 → 注销biocrf（！！！Warnings，dangerous to unix system management！！！）

# 用户被注销，与AVITI机器连接的原smb协议断开
# 添加新的samba用户
sudo smbpasswd -a biocrfgz
# 修改密码为biocrfgz
```shell