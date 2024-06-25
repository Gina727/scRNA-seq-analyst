
<h1 align="center">
  <br>
  <img src="https://github.com/Gina727/scRNA-seq-analyst/blob/main/app-api/scRNA%20logo.png" width="400">
  <br>
   scRNA-seq Analysis Copilot
  <br>
</h1>

<h4 align="center">A single cell sequencing data analysis interactive helper. Learn scRNA-seq Analysis through natural language interaction with AI charbot.</h4>

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#credits">Credits</a> •
  <a href="#related">Related</a> •
  <a href="#license">License</a>
</p>

## Key Features

* Answering questions with knowledge
  - 9 preloaded papers with 658,018 words were uploaded for retrieval and answer questions related to general ideas, workflow, sequencing technology, annotation, and mitochondrial genes.
* Single cell sequencing analysis workflow
  - The copilot is capable of loading datasets from given direct url links, process them with API R functions on the server end.
* Graph demonstration
  - The copilot shows generated plots to the users directly.

## How To Use
This application is supported by and published on Dify platform. Users might access the agent via <a href="http://dify.docai.net/chat/kyCD408hEc5p17yu">Ready-to-use AI WebApp</a>.
To clone and run this application, you'll need [Dify](http://dify.docai.net/apps) for GPT agent and [R](https://www.r-project.org/) for writing R code and build API installed on your server.

Part One. Build your server
1. Connect to your server through SSH applications like PuTTY.
2. Enter command lines on your server:
```bash
# Clone this repository
$ git clone https://github.com/Gina727/scRNA-seq-analyst

# Go into the repository
$ cd scRNA-seq-analyst
```
2. Build the image of app-api
```bash
$ docker build -t new-api .
$ docker compose up -d
```
You may use $ docker ps to check for running containers.

3. Add folders, rds and im
4. Open "http://your-server's-address:8000/openapi.json"
5. Copy the content of the webpage.

Part Two. Build your agent
1. Register on Dify
2. Create an agent
3. Enter the prompt from this repository: https://github.com/Gina727/scRNA-seq-analyst/blob/main/prompt.txt
4. Load papers in .pdf file from the knowledge list: https://github.com/Gina727/scRNA-seq-analyst/blob/main/knowledge_list.txt to the "Knowledge" section.
5. Go to Tools section. Create new API tools by pasting the openapi.json in the corresponding section.
6. Done!

You may edit the "new_api.R" R script to add new functions.

## Credits

This software uses the following open source packages:

- [Plumber](https://github.com/rstudio/plumber)
- [Seurat](https://github.com/satijalab/seurat)
- [Seurat Data](https://github.com/satijalab/seurat-data)
- [ScType](https://github.com/IanevskiAleksandr/sc-type)

## License

MIT

---

