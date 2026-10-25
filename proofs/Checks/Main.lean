module

import Lean
import Checks.Examples

open Lean

namespace ProofAudit

/-- Export structural expressions, not human pretty-printing that could hide a changed type. -/
def declarationJson (env : Environment) (sourceModule name : Name)
    (info : ConstantInfo) : IO Json := do
  let axioms ← (collectAxioms name : CoreM (Array Name)).toIO'
    { fileName := "proof-audit", fileMap := default } { env := env }
  let doc ← findDocString? env name
  -- Direct ranges do not borrow a parent's location for generated recursors or
  -- auxiliary proofs. Structure constructors still have direct ranges, so identify
  -- their generated construction rule using Lean's structure metadata as well.
  let ranges ← (findDeclarationRangesCore? name : CoreM _).toIO'
    { fileName := "proof-audit", fileMap := default } { env := env }
  let structureConstructor := match info with
    | .ctorInfo ctor => isStructure env ctor.induct
    | _ => false
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
    ("authored", toJson (ranges.isSome && !structureConstructor)),
    ("kind", toJson kind),
    ("type", toJson (reprStr info.type)),
    ("value", toJson (info.value? (allowOpaque := true) |>.map reprStr)),
    ("uses", toJson (info.getUsedConstantsAsSet.toArray.map Name.toString)),
    ("axioms", toJson (axioms.map Name.toString)),
    ("doc", toJson doc)]

/-- Load the actual built mathematical modules and export their checked declarations. -/
unsafe def metadata : IO Json := do
  initSearchPath (← findSysroot)
  enableInitializersExecution
  let inventory ← IO.FS.readFile "theorems.json"
  let json ← IO.ofExcept (Json.parse inventory)
  let modules ← IO.ofExcept (json.getObjValAs? (Array String) "modules")
  let imports : Array Import := modules.map fun source =>
    { module := (source.dropEnd 5 |>.toString.replace "/" ".").toName, importAll := true }
  let env ← importModules imports {} (loadExts := true)
  let mut declarations := #[]
  for (name, info) in env.constants.toList do
    if let some idx := env.getModuleIdxFor? name then
      let sourceModule := env.header.moduleNames[idx]!
      -- Ownership follows every maintained module, including checker modules that a
      -- model might import. Keep generated constants for transitive fingerprints too.
      if imports.any (·.module == sourceModule) then
        declarations := declarations.push (← declarationJson env sourceModule name info)
  return Json.mkObj [
    ("lean_version", toJson Lean.versionString),
    ("examples", NetGraphCore.Checks.caseJson),
    ("declarations", toJson declarations)]

end ProofAudit

/-- Run the packaged audit with metadata produced from Lean's real environment.
Python supplies SHA-256, filesystem and JSON checks; it does not prove theorems. -/
public unsafe def main (args : List String) : IO UInt32 := do
  let data ← ProofAudit.metadata
  if args == ["--metadata"] then
    (← IO.getStdout).putStrLn data.compress
    return 0
  let python := (← IO.getEnv "PROOFS_PYTHON").getD "python3"
  let out ← IO.Process.output {
    cmd := python
    args := #["-I", "-B", "Checks/audit.py"] ++ args.toArray
  } (some data.compress)
  (← IO.getStdout).putStr out.stdout
  (← IO.getStderr).putStr out.stderr
  return out.exitCode
