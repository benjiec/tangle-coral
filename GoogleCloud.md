## Google Cloud Setup

Create project

Enable services
  * Batch API
  * Artifact Repository API
  * Compute Engine API
  * BigQuery API

Create service account for nextflow, e.g. nextflow-service-account@needle-489321.iam.gserviceaccount.com

Add roles to service account
  * Artifact Registry Reader
  * Batch Agent Reporter
  * BigQuery Data Editor
  * Storage Object Admin

Also, create a key for this service account, download the JSON key and put in ~/.config/gcloud

Some initialization stuff

```
gcloud init
gcloud auth application-default login
gcloud auth configure-docker us-east1-docker.pkg.dev
```

Create a docker repository

```
gcloud artifacts repositories create needle-docker \
    --repository-format=docker \
    --location=us-east1 \
    --description="Docker repository for Needle Nextflow workflow"
```

Create a bucket for storing files

```
gcloud storage buckets create gs://needle-files\
    --location=us-east1 \
    --uniform-bucket-level-access
```

And copy some files

```
gcloud storage cp -r ./kegg-downloads gs://needle-files/kegg-downloads
gcloud storage cp -r data/ko_thresholds.tsv gs://needle-files/kegg-downloads/
gcloud storage cp -r ./pfam-downloads gs://needle-files/pfam-downloads
```


## Building Images

```
docker build --platform linux/amd64 -t us-east1-docker.pkg.dev/needle-489321/needle-docker/needle-core:latest .
```

We can create a repository, then push the image

```
docker push us-east1-docker.pkg.dev/needle-489321/needle-docker/needle-core:latest
```

Test the image

```
gcloud batch jobs submit test-job-v1 \
    --location us-east1 \
    --config - <<EOF
{
  "taskGroups": [
    {
      "taskSpec": {
        "runnables": [
          {
            "container": {
              "imageUri": "us-east1-docker.pkg.dev/needle-489321/needle-docker/needle-core:latest",
              "commands": ["python3", "scripts/classify/classify.py"]
            }
          }
        ]
      },
      "taskCount": 1
    }
  ],
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}
EOF
```

Then you can look at the status of the job

```
gcloud batch jobs describe test-job-v1 --location us-east1
gcloud batch jobs delete test-job-v1 --location us-east1
```
