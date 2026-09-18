# Plan Mode Rules

- When generating or saving plans in Plan Mode, you MUST write the plan file directly into the project directory under `./.claude/plans/`.
- Do NOT save plans to the global default `~/.claude/plans/` directory.
- When instructed to execute a plan from  `./.claude/plans/` directory, use that plan only. Ignore all other plans in `./.claude/plans/`.
- Optimize the plan instructions for  prompting, especially caching 
  efficiency - group the needed files and prepend `@` to each file, stripping down the quotes if present,
  for example (this example is not from this repo):
```txt
  > I'm refactoring the payment processing. The current flow is in:
  @src/services/payment.ts (main logic)
  @src/api/stripe.ts (payment provider integration)
  @src/models/Transaction.ts (data model)
```
