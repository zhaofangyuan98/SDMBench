
# Docker

This Docker repository is continuously updated as more new spatial clustering methods are published in the field.

## 1. Constructing a New Docker Image and Publishing to Docker Hub

### 1.1 Building a New Docker Image on the Local Machine

The images we provide are based on nvidia/cuda:11.7.1 and ubuntu18.04. If you need to change the cuda and system versions, you can modify them by modifying the dockerfile we provide. For users who wish to utilize the image provided directly, please proceed to step 2. For those who wish to construct a custom Docker image, please ensure a Dockerfile exists within the directory of interest. The following command can then be executed to build the Docker image:

```bash
sudo docker build -t <custom_image_name> .
```

This command will leverage the Dockerfile to construct a Docker image named my_custom_image.

### 1.2 Publishing the Docker Image to Docker Hub or Other Docker Image Repositories

Users should initially log into Docker Hub or their preferred Docker image repository on their local machine. For instance, to log into Docker Hub, the following command can be used:

```bash
sudo docker login
```

It is required for users to create a Docker Hub account in advance and authenticate using their respective username and password.

### 1.3 Publish the Docker Image to Docker Hub

Users are advised to substitute <username> with their Docker Hub username in the following commands:

```bash
sudo docker tag <custom_image_name>:latest <username>/<custom_image_name>:latest
sudo docker push <username>/<custom_image_name>:latest
```

## 2. Deploying a Docker Image on a New Server

### 2.1 Install Docker on the New Server

The process of installing Docker may differ based on the server's operating system. The following commands illustrate the basic steps to install Docker on Ubuntu:

```bash
sudo apt-get update
sudo apt-get install docker-ce
```

### 2.2 Pull the Image from Docker Hub
  
| Username/Image_name | Tag |
|:-------:|:-------:|
| sdmbench/leiden | latest |
| sdmbench/louvain | latest |
| sdmbench/stagate_gpu | latest |
| sdmbench/spaceflow_gpu | latest |

We provide the image above：

eg. sudo docker pull linsenlin/spatial_clustering_gpu:latest

```bash
sudo docker pull <username>/<custom_image_name>:latest
```

## 3. Launch the Docker Container on the New Server

At this point, users can initiate their Docker container on the new server. The following command will run the container and map port 8888 of the server to port 8888 of the container:

eg. sudo docker run -p 8888:8888 sdmbench/stagate_gpu:latest

```bash
sudo docker run -p 8888:8888 <username>/<custom_image_name>:latest
```

If user's device configuration is based on nvidia/cuda:11.7.1 and ubuntu18.04, and the user wishes to utilize GPU, they should execute the following command:

eg. sudo docker run --gpus all -p 8888:8888 sdmbench/stagate_gpu:latest  

```bash
sudo docker run --gpus all -p 8888:8888 <username>/<custom_image_name>:latest
```

## 4. Access Jupyter Notebook

Upon running the Docker container on the new server, users should be able to access Jupyter Notebook via port 8888 on the server. If the server's IP address is 'server_ip', users can enter the following URL into their web browser:

`http://server_ip:8888`
 
This guide provides users with the necessary steps to construct, publish, and deploy Docker images to enhance their software deployment process.
