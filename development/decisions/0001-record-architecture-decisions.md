# 0001. Record Architecture Decisions

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user and Claude. This follows the convention shared with the laboratory and the sibling libraries.

## 1. Context

calx is long-lived, and other libraries depend on it. Later sessions need to know why its interfaces are the way they are, including sessions that do not remember the original discussion.

## 2. Decision

- Record each decision that shapes calx's structure, interfaces or tooling as a numbered record in `development/decisions/`.
- Never rewrite a record. A changed decision gets a new record that supersedes the old one.

## 3. Alternatives Considered

- **Commit messages only.** Rejected: the reasons are hard to find later.
- **A single design document, edited in place.** Rejected: it loses the history of why things changed.

## 4. Consequences

- Each decision costs a small amount of writing.
- `docs/architecture.md` describes the current state. The records explain how that state came about.
