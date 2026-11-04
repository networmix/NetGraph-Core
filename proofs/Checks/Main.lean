module

import Lean
import Checks.Examples

open Lean

namespace ProofAudit

/-- Export declaration ownership and actual transitive axioms from the built environment. -/
def declarationJson (env : Environment) (sourceModule name : Name)
    (info : ConstantInfo) : IO Json := do
  let axioms ← (collectAxioms name : CoreM (Array Name)).toIO'
    { fileName := "proof-audit", fileMap := default } { env := env }
  let kind := match info with
    | .thmInfo _ => "theorem"
    | .defnInfo _ => "definition"
    | .axiomInfo _ => "axiom"
    | .opaqueInfo _ => "opaque"
    | .inductInfo _ => "inductive"
    | .ctorInfo _ => "constructor"
    | .recInfo _ => "recursor"
    | .quotInfo _ => "quotient"
  return Json.mkObj [
    ("name", toJson name.toString),
    ("user_name", toJson (privateToUserName name).toString),
    ("source_module", toJson sourceModule.toString),
    ("kind", toJson kind),
    ("axioms", toJson (axioms.map Name.toString))]

/-- Load the actual built mathematical modules and export their checked declarations. -/
unsafe def metadata : IO Json := do
  initSearchPath (← findSysroot)
  enableInitializersExecution
  let python := (← IO.getEnv "PROOFS_PYTHON").getD "python3"
  let listed ← IO.Process.output {
    cmd := python
    args := #["-I", "-B", "Checks/audit.py", "--modules"]
  }
  if listed.exitCode != 0 then
    throw (IO.userError listed.stderr)
  let json ← IO.ofExcept (Json.parse listed.stdout)
  let modules ← IO.ofExcept (fromJson? json : Except String (Array String))
  let imports : Array Import := modules.map fun name =>
    { module := name.toName, importAll := true }
  let env ← importModules imports {} (loadExts := true)
  let mut declarations := #[]
  for (name, info) in env.constants.toList do
    if let some idx := env.getModuleIdxFor? name then
      let sourceModule := env.header.moduleNames[idx]!
      -- Include every local module, including checker modules a model might import.
      -- Generated declarations also need their transitive axioms checked.
      if imports.any (·.module == sourceModule) then
        declarations := declarations.push (← declarationJson env sourceModule name info)
  return Json.mkObj [
    ("lean_version", toJson Lean.versionString),
    ("examples", NetGraphCore.Checks.caseJson),
    ("declarations", toJson declarations)]

end ProofAudit

/-- Run the packaged audit with metadata produced from Lean's real environment.
Python checks source safety, axioms and exact examples; it does not prove theorems. -/
public unsafe def main (args : List String) : IO UInt32 := do
  let data ← ProofAudit.metadata
  let python := (← IO.getEnv "PROOFS_PYTHON").getD "python3"
  let out ← IO.Process.output {
    cmd := python
    args := #["-I", "-B", "Checks/audit.py"] ++ args.toArray
  } (some data.compress)
  (← IO.getStdout).putStr out.stdout
  (← IO.getStderr).putStr out.stderr
  return out.exitCode
