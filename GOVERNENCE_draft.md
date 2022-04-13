<!-- This document is borrowed from ggplot2 governance doc, which was:  heavily adapted version of
the Benevolent dictator governance model by Ross
Gardler and Gabriel Hanganu licensed under a Creative Commons
Attribution-ShareAlike 4.0 International License. -->

# Omicser Contribution and Governance Policies

This document describes the contribution process and governance policies of the NDCN Omicser project.  As this project aims to solve for a broad set of needs across the NDCN community there will need to be honest and thoughtful discussion around features and priorities.

## Contribution Process
Before making a contribution, please take the following steps:
1. Check whether there's already an open issue related to your proposed contribution. If there is, join the discussion and propose your contribution there.
2. If there isn't already a relevant issue, create one, describing your contribution and the problem you're trying to solve.
3. Respond to any questions or suggestions raised in the issue by other developers.
4. Fork the project repository and prepare your proposed contribution.
5. Submit a pull request.

NOTE: All contributors are implicitly sharing open license to the code submitted.  See also our [CONTRIBUTION document](https://github.com/ndcn/omicser/blob/main/CONTRIBUTING.md)


## Governance

### Roles and responsibilities

This project has a large community of ***users*** and ***contributors***, a team of ***maintainers***, the ***”Brain Trust”***, and a ***NDCN Github Org Lead***.

#### Users

People who browse and visualize -omics data with the _NDCN Omics Browser_ (or `omicser`) are the most important members of the community; without these users, this project would have no purpose.

Users are encouraged to participate in the life of the project and the community as much as possible. User contributions help ensure that the project is satisfying users' needs. Common user activities include (but are not limited to):

- evangelizing about the project
- asking and answering on community forums
- providing moral support (a 'thank you' goes a long way)

Users who continue to engage with the project and its community will often find themselves becoming more and more involved. Such users may then go on to become contributors, as described above.

#### Contributors

Contributors interact with the project on GitHub by filing new issues, improving existing issues, or submitting pull requests. Anyone can become a contributor: there is no expectation of commitment to the project, no required set of skills, and no selection process. The only obligation is to follow the [code of conduct](CODE_OF_CONDUCT.md).

Specific advice for contributing to the project can be found in [CONTRIBUTING.md](CONTRIBUTING.md).

#### Maintainers

Maintainers are collectively responsible for day-to-day development of the package, including responding to issues and reviewing pull requests. They are GitHub administrators and [package authors](https://github.com/ndcn/omicser/blob/master/DESCRIPTION#L5), which means that they have the ability to make changes to project code, and receive credit when others cite the package.

While maintainers can modify code directly, this ability is rarely used. Instead, changes are proposed as pull requests, and are only merged after they have been reviewed by at least one other core developer. Changes to the API (especially breaking changes) must also be approved by the project lead.

Maintainers are recruited from contributors. An invitation to join the core team can be extended to anyone who has made a major contribution, either through a small number of large changes, or a consistent pattern of smaller contributions. Any existing core developer can propose a contributor be invited to the core team by contacting the project lead. The project lead will the confirm the invitation with the other maintainers.

The maintainers of omicser are:

* [Rico Derks](https://github.com/rderks) Head Maintainer
* [Andy Henrie](https://github.com/ergonyc) Original Author

Our _Head_ Maintainer, [Rico Derks](https://github.com/rderks) will additionally delegate tasks, communicate proactively with the project lead, and ensure that pull requests, issues, etc are handled in a timely manner.

All maintainers are bound by the [code of conduct](CODE_OF_CONDUCT.md).

<!-- More details can be found in the [maintainers guidelines](MAINTAINER_GUIDELINES.md).-->

#### “Brain Trust”
The NDCN Browser  ***“Brain Trust”*** is the group of NDCN  PIs, their lab members, and CZI personnel who continue to steer the open source project.  The “Brain Trust” is responsible for:

* Specify and evolve the strategic objectives of the project.
* Enabling the project survives in the long term.

Together the group of individuals who constitute the ***“Brain Trust”*** will proactively articulate the overall strategies, goals and objectives of the project.  Although the ***”Brain Trust*** is everyone involved in the bi-monthly discussions, the primary responsibility is shared among the PIs of that group.   Currently, all of the ***Maintainers*** are __also__ memebers of the ***Brain Trust***.  They will ensure that the ***(Head) Maintainers*** tactical objectives are aligned with the overall project.  

#### Project NDCN Github Organization lead
The ***NDCN Github Organization Lead***  leads,  [Chris Sifuentes](https://github.com/cjsifuen), is responsible for:
* Clearly communication of the NDCN Github Organization policies and objectiveson the strategic objectives of the project.
* Mediating any conflicts amongst the maintainers in repositories within the NDCN Github Organization.

The project lead is bound by the [code of conduct](CODE_OF_CONDUCT.md).

### Decision-making process

This project makes decisions according to a consensus model where suggestions are considered and discussed between the community and maintainers, typically in GitHub issues. If the community questions a decision, the ***“Brain Trust”*** and/or the ***NDCN Github Organization Lead*** may review it and either uphold or reverse it.
