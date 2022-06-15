# Git-based dependency management workflow #

## Introduction ##
We present a **git-based dependency management workflow**. It is composed of:
1. a way to organize git repositories for developing multiple projects sharing common libraries
3. a set of bash, git and cmake tools to help implementing the workflow.

The main driving factor for using git-based dependency management comes from the fact that we want to have small-sized semi-independent libraries while being able to develop in several of them seamlessly (modifying files here and there without worrying too much if it comes from one library or another). Since all libraries are being developed at the same time, we need to keep track of library versions. But we don't want to keep track of the versions in a heavy-weight fashion, meaning we don't want to manually specify a new version number for a library each time we do a git commit (as we would do with e.g. SemVer): on a day-to-day basis, identifying library versions by their commit number should be enought and not need any additional step.

To give a more concrete idea of inter-dependent libraries we want to tackle, see the example of dependencies between libraries on the graph below.

![Dependency graph of connected git repositories](./project_dep_graph.svg)


## Scope ##

The projects we are talking about are **fine-grained**: a few people working on a set of small libraries. We picture 2/3 people working on the same library at the same time (if the number of people is greater, maybe it means that the library is too big and should be split into several ones).

We think SemVer-based version management is also useful and should also be used, but on a more coarse-grained time frame, for sets of bigger libraries, or for cross-team development.


## Features and goals ##
1. Versioning done by git. Git has become the de-facto standard and has gained a large consensus as being the right tool for the job.
2. Ability to **easily create a library** as we see fit, or extract a factored piece of code into a library. A well-written library has a well-defined role, a weak coupling with its dependencies, and no coupling at all to the other projects that depend on it. Thus, it is easier to learn, share, maintain and test.
3. **Each library has its own git repository** and can be built by itself. This simplifies testing. Furthermore, it also simplifies learning it and ensures there is no hidden coupling with a "container" repository.
4. Possibility to develop collaboratively on several projects and library repositories at the same time.
5. Dependencies between repositories (both library and applications) are **explicit, but lightweight**. This means that when working on several repositories, we do not manually create version numbers for each dependency repository and then tell the dependent repository which version number it should use. This would be perfectly fine at a higher level of granularity, but here we are at a very fine level and versionning should be commit-based: "libA at commit 5f25414f works with libB at commit f16af47b", that's all: we do not ensure any further compatiblity.
6. Possibility to depend on libraries not integrated within this framework of interdependent git repositories. Possibility to add other version-tracking mechanisms (e.g. SemVer).
7. No forced updates. While merging with the latest developments should be done as soon as it is possible, staying on the same commit of a dependency and not updating it while doing something else is perfectly acceptable.

### Limitations ###
1. API stability. While it is certainly a good idea to ensure API stability and proper versioning for top-level projects, there is a tradeoff with flexibility. Here, we do not want to garantee API stability regarding small-scale libraries as this would hamper the ability to easily split projects into libraries for the sake of ensuring API stability for each new library.
2. Complexity. Interdependent libraries developed collaboratively on different versions and git branches is a complex topic. There is no magic here.
3. User-friendliness. Git submodules come naked, and while we try to provide some convenience scripts, tools are not always complete or well-polished.


## Description ##
### Folder structure ###

The project structure of `My_project` given in the dependency graph above would look like this:

```
My_project/
├── README.md
├── LICENSE.md
├── CMakeLists.txt        # Top-level CMakeLists.txt for compiled project.
│
│
├── My_project/           # Code of the library. No include/src separation (not useful). Contains unit tests.
├── external/             # Library dependencies. They are submodules.
│   ├── Maia/
│   ├── Tasky/
│   ├── std_e/
│   └── project_utils/
│
├── doc/
├── examples/
│
├── scripts/              # Additional scripts (building, building tests,
│                         #                deployment machine-specific instructions...)
├── test/                 # Functionnality/acceptance/performance/deployment tests
└── ...                   # other sub-folders
```

The important part is that all the libraries that `My_project` depends upon are git submodules placed in the sub-folder `external`.

### Dependency management ###
There are three kinds of dependencies:
* **Project-related** dependencies
* **Other-project**  dependencies (as per feature 6.)
* **System** dependencies

The latter two are taken care of through CMake (or other build system tools if not possible) and proper versionning. Regarding project-related dependencies, they must be managed at the commit level but we are not aware of any git-based (or other version-control system) dependency managers. Our compromise is to use git submodules with a disciplined approach. In our example of library dependencies, dependencies form a direct acyclic graph of depth 4. It is complicated:
1. It is 4 levels deep, which mean that dependencies (example: `Project_B` depends on `Maia`) have themselves dependencies (`Maia` depends on `std_e`)
2. It is NOT a tree. For `My_project` depends on `Maia` and `Tasky`, and both depend on `std_e` ("diamond-shaped dependency"). Which raises two questions:
    1. What should we do if `Maia` and `Tasky` depend on two different versions of `std_e`?
    2. Should we include the content of `std_e` twice? If yes, when changing it during the development of `My_project`, which one do we change?

#### Chosen solution ####
1. Depending on two versions of the same library is forbidden/not supported.
2. Each repository `r` must have all its project-related dependencies `d` taken into account as git submodules placed in `root_of_r/external/d`.
3. Each repository `r` **also must have its indirect dependencies `i_d` stored in `root_of_r/external/i_d`**.
4. When working on repository `r`, *only its submodules are checked out*. The submodules of the submodules (indirect dependencies) are already checked out at the first level, and they are supposed to have the same version. Hence there is no need to have them twice.

* [+] The approach is relatively simple because from a working tree perspective, there is only one submodule depth.
* [-] The top-level repository must know all its indirect dependencies. If we were to use a complex layering of many libraries, this could be a problem, but:
    * In reality, the total number of repositories involved is small. If it grows too much, then some *project-related* dependencies should be separated as *other-project* dependencies.
    * In our example, indirect dependencies are also direct dependencies, so it doesn't change anything. Our experience is that we are in this scenario most of the time.
* [+] At this fine-grained "white-box" level, it would not make sense to depend on two versions of the same library. If this is the case at some point, it happened because of an update at one side. So either the discrepancy should be fixed, or the update should be reverted.
* [-] We must configure local git repositories in order to see new commits of `std_e` for all its dependencies (see section "Git additionnal commands")

## Workflow ##
* Put the contents of `scripts/git/submodule_utils.sh` in your **bash profile**. It will configure git to print more info if a change is made to a submodule. The git aliases `sclone`, `scheckout`, `spull` and `spush` will also be available.
* **Clone a repostory** with `git sclone <repo-name>`. This will also make sure that submodules required several times in the dependency graph are only cloned once and share the same repository.
* When you **switch to a branch**, use `git scheckout` so that the dependencies' working trees will be updated to the versions referenced by the current branch of main repository.
* When you are **pulling**, use `git spull` to update the dependencies. If you pull with `git pull` or `git fetch + git merge`, the dependencies' working trees are still on older versions. To update them to the versions referenced by the current branch of main repository, use `git submodule update`.
* If you develop only in the main project, then you don't need to care about submodules for anything else.
* If you made a **change to a submodule** `std_e` (located in `<main-repo-path>/external/std_e`) then
    * You should see it with e.g. `cd <main-repo-path>; git status`.
    * You can commit the change in the submodule by going into its folder (`cd <main-repo-path>/external/std_e`)
    * When in the submodule folder, git does as if it were a regular git repository.
    * Be aware that, by default, submodules are on a "detached head" state (see `git status`) which mean you are not on a branch (not even master).
        * If you are to make commits, you should go on a branch.
        * You can ask git for the branch of the current commit with `git rev-parse HEAD | xargs git name-rev`.
        * Say the branch is `master`, you can go on it by `git checkout master`.
    * You can then commit your change (`git add ...; git commit ...`).
    * If `std_e` is also a dependency of another submodule `Maia`, then changes to `std_e` will also be reported as changes to `Maia`. After having commited the changes to the `std_e`, you can commit the update of `std_e` into `Maia`, and then the update of `std_e` and `Maia` into the main repository.
* When you are **pushing** a repository, if you changed one of its submodules (new commits), then make sure to also push it (`git spush` will do it for you).
* TODO create git alias to automate commits of just submodule dependency updates

## Common cases ##

The example repository names still refer to the dependency graph above.

### Updating submodules ###
* Say that we are working on `My_project`. If submodule `std_e` has been changed outside of `My_project` (e.g. through developpers working on `Project_B`):
    * **Most of the time, it doesn't matter**. Don't do anything special regarding `std_e`. Use `git spull` on `My_project` to get a coherent, new version of `My_project`. It will **not** pull the latest `std_e` changes created by the unrelated `Project_B`. This is the correct behavior, because it ensures that a particular commit of `My_project` is not silently affected by new versions of its dependencies (here, by a new version of `std_e` developped in a different context than `My_project`).
    * If you want to update `std_e`, go to `external/std_e` and pull. Then when you come back to the main folder project, you should see that `std_e` has an updated version. Commit the change. You may have to commit changes to other submodules that depend on `std_e` (e.g. `Maia`) before that.

### Modifying a submodule ###
* Say that we are working on `My_project` and in the development process, we want to develop a functionality that is more in the scope of `std_e`:
    * We do the development that modifies the source files of `std_e`.
    * Then `cd My_project/external/std_e && git commit && git push`.
    * Then `cd ../.. && git status`. We should see that `std_e` now has a new version.
    * We can update `My_project` to reflect on the fact that we are now using a new version of `std_e`. It is done with the standard `git commit` command. From the point of view of the `My_project` git repository, the only thing that is modified is the commit number of `std_e`

**Do not forget to push** the development in `std_e` when you push the development in `My_project`. If you forget, a user may update `My_project`, and find that the `std_e` commit that it references is not known by the server (because you did not push it).

### Switching branches - managing DETACHED_HEAD ###
* Use `git scheckout` instead of `git checkout` to switch branches. If you don't, you will switch from the old to the new branch, but the dependencies will stay at the versions that were used by the old branch.
* When you switch `My_project` to a new branch `my_branch` with `git scheckout`, `std_e` will be placed at the commit number that was registered by `My_project/my_branch`. However, that commit number could be anywhere in the git history, e.g. `master~10` (git convention that means "10 commits before the current master"). Since `master~10` is not the tip of the branch, git places `std_e` in a `DETACHED_HEAD` state.
* A `DETACHED_HEAD` state is a way for git to say: "if you make a commit, I don't know on which branch to put it". For example, if we are at `master~10`, what is the branch of the next commit? It can't be `master` since it already has a next commit, `master~9`!
* One solution is then to create a branch from `master~10`: `git checkout -b new_std_e_branch`
* Another solution is to go to the tip of the branch: `git checkout master`. Be aware that you may need to update the sources of `My_project` since you changed to a new version of its dependency `std_e`.

Note: actually `git status` does not tell you that you are on commit `master~10`, it just gives the commit hash. You can find the branch by looking the history with a GUI or you can use `git rev-parse HEAD | xargs git name-rev`.

### Adding submodules ###
* Let's say you are creating the library `Tasky` and want to use `std_e` as a submodule. Then, as explained in (see [this section](#chosen-solution)), you **must** also add `project_utils` as a submodule.
* Notice however that when adding `std_e`, `project_utils` will be installed as a submodule of `std_e`. You then have `project_utils` sources two times: in `Tasky/external/project_utils`, and in `Tasky/external/std_e/external/project_utils`. You want to delete the content of the second directory to replace it by a git symbolic link to the first directory. For that:
    * `cd Path/To/Tasky/external/std_e/external/project_utils && rm -rf *`
    * `echo "gitdir: Path/To/Tasky/.git/modules/external/project_utils" > .git`


## Provided commands and Git aliases ##
### sclone ###
`git sclone` should be used instead of `git clone` when working with this framework. It will download and configure submodules as expected by the framework.

#### details ####
We need `git sclone` in order to handle "diamond-shaped dependencies". For example, suppose the following scenario:
 * working on `My_project` (at version v0), with dependencies `Maia` (at v0) and `std_e` (at v0).
 * note that even if both `My_project` and `Maia` depend on `std_e`, there is only one `std_e` repository, located at `My_project/external/std_e` (see [this section](#chosen-solution))
 * we modify `std_e`, submodule of `My_project`, then commit the change. The git working tree of `std_e` is now at v1.

Now for both `Maia` and `My_project`, the dependency `std_e` needs to be seen has changed to v1. For `My project`, this works out-of-the box since `std_e` is located at `My_project/external/std_e`. But for Maia, we would expect a kind of symbolic link from `My_project/external/Maia/external/std_e` to `My_project/external/std_e` to be informed, when inside `My_project/external/Maia`, that `std_e` has been modified. This is actually exactly what is put in place by `git sclone` (the git equivalent of a symbolic link actually being that `My_project/external/Maia/external/std_e/.git` contains the line `gitdir: My_project/.git/modules/std_e`).

### git\_config\_submodules ###
`git sclone` is in reality mostly a chaining of:
* `git clone`
* `cd $cloned_repo`
* `git submodule update --init`
* `git_config_submodules`
