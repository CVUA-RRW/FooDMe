# Contribution guide

Maintaining a repository is a lot of work and it is only possible with the
input for its users. That's why we love to get your feedback and continuously
improve! We want to contibuting to this project as easy and transparent as
possible. That's why we ask all contributions to go through a GitHub workflow.

## Ask a question, report a bug, or suggest new features

Wether you wish to report a bug, discuss the current state of the code,
how to use it, or propose new features, the GitHub repository is the place to
start.

Go to 'Issues' in the repository's menu and use the search bar to look for your
issue/question, maybe the discussion already exists! if you don't find what
you are looking for, use the green button
['New Issue'](https://github.com/CVUA-RRW/FooDMe/issues)
and select the correct template (Question, Bug report or New feature).

### Bug reports

Solving Bugs is not always easy and usually requires to be able to precisely
understand what happened. Try to include the following in your report, so we
can solve the problems quickly:

- A quick summary and/or background
- The software versions
- A precise description of the problem
- Join input data that reproduce the problem if possible
- Attach relevant log files to the issue

## Improving the documentation

Writting a documentation is a big task and we are gratefull for any help
to expand and improve it.

If you find a specific section of the documentaiton unclear, or would like to
see it in other languages, we would love to hear about it, or if you feel like
it, try suggesting modifications.


## Writting tests

The critical parts of the workflow are currently verified by unit tests.
These could always be improved and/or expanded. If you like to write tests,
check the `.tests/unit` folder and suggest changes.

## Submitting changes

All changes in the repository are made through pull requests.

### Github worflow

To submitt a pull request:

- Fork the repository and create your branch from `main`
- Make your changes to the code base/documentation/tests
- If you have added code that shoul dbe testes, add tests
- If you have changed the parameters or formating of the inputs, update the documentation
- Ensure that your code passes the tests and lints, this will be automatically
  tested when you submit your pull request
- Issue your pull request and wait for a maintainer to review the changes

Please provide precise and sufficient information on the changes you
performed when you submit your pull request.

### Tests

It is a good idea to locally run tests before submitting your pull request.

#### Unit tests

We use pytest for running unit tests. Make sure your current environment supports
a recent version of python (3.9 or 3.10).

```bash
python -m pip install --upgrade pip
pip install pytest
pip install -r .tests/unit/requirements.txt  # Install tests dependencies
pytest .tests/unit
```

#### Integration test

It is always a good idea to verify that the workflow runs properly as a whole.
For this activate a conda enviroment with a recent version of snakemake (see the
documentation) and run the following test:

```bash
snakemake --cores 1 --use-conda --configfile .tests/integration/config/config.yaml
```

### License

By submitting changes to this repository, you agree that your contributions
will be licensed under its [BSD-3 Clauses License](LICENSE).
