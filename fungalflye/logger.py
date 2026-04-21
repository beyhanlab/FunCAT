"""
funcat.logger
~~~~~~~~~~~~~~~~~
Comprehensive logging system for FunCAT runs.
Captures user selections, commands, timing, and results.
"""

import json
import time
from pathlib import Path
from datetime import datetime


class FuncatLogger:
    """Comprehensive logger for FunCAT pipeline runs."""
    
    def __init__(self, outdir):
        self.outdir = Path(outdir)
        self.log_file = self.outdir / "funcat_run.log"
        self.start_time = time.time()
        self.run_data = {
            "start_time": datetime.now().isoformat(),
            "version": "0.3.0",
            "user_selections": {},
            "commands_executed": [],
            "module_results": {},
            "timing": {},
            "errors": [],
            "final_stats": {}
        }
        
        # Ensure output directory exists
        self.outdir.mkdir(exist_ok=True)
        
        self._write_header()
    
    def _write_header(self):
        """Write initial log header."""
        header = f"""
{'='*80}
FunCAT Pipeline Run Log
Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
{'='*80}

"""
        with open(self.log_file, 'w') as f:
            f.write(header)
    
    def log_user_selection(self, step, selection, value):
        """Log user input selections."""
        self.run_data["user_selections"][step] = {
            "selection": selection,
            "value": value,
            "timestamp": datetime.now().isoformat()
        }
        self._append_log(f"USER INPUT: {step} = {selection}: {value}")
    
    def log_command(self, command, cwd=None):
        """Log command execution."""
        timestamp = datetime.now().isoformat()
        cmd_entry = {
            "command": command,
            "timestamp": timestamp,
            "cwd": str(cwd) if cwd else None
        }
        self.run_data["commands_executed"].append(cmd_entry)
        
        cwd_str = f" (cwd: {cwd})" if cwd else ""
        self._append_log(f"COMMAND: {command}{cwd_str}")
    
    def log_module_start(self, module_name):
        """Log start of a pipeline module."""
        timestamp = datetime.now().isoformat()
        self.run_data["timing"][f"{module_name}_start"] = timestamp
        self._append_log(f"MODULE START: {module_name}")
    
    def log_module_end(self, module_name, result_path=None, stats=None):
        """Log end of a pipeline module."""
        timestamp = datetime.now().isoformat()
        self.run_data["timing"][f"{module_name}_end"] = timestamp
        
        # Calculate module duration
        start_key = f"{module_name}_start"
        if start_key in self.run_data["timing"]:
            start_time = datetime.fromisoformat(self.run_data["timing"][start_key])
            end_time = datetime.fromisoformat(timestamp)
            duration = (end_time - start_time).total_seconds()
            self.run_data["timing"][f"{module_name}_duration"] = duration
        
        # Store module results
        self.run_data["module_results"][module_name] = {
            "result_path": str(result_path) if result_path else None,
            "stats": stats,
            "completed_at": timestamp
        }
        
        duration_str = f" ({self.run_data['timing'].get(f'{module_name}_duration', 0):.1f}s)" if f"{module_name}_duration" in self.run_data["timing"] else ""
        result_str = f" → {result_path}" if result_path else ""
        self._append_log(f"MODULE END: {module_name}{duration_str}{result_str}")
        
        if stats:
            for key, value in stats.items():
                self._append_log(f"  {key}: {value}")
    
    def log_error(self, error_msg, exception=None):
        """Log errors and exceptions."""
        timestamp = datetime.now().isoformat()
        error_entry = {
            "message": error_msg,
            "exception": str(exception) if exception else None,
            "timestamp": timestamp
        }
        self.run_data["errors"].append(error_entry)
        
        self._append_log(f"ERROR: {error_msg}")
        if exception:
            self._append_log(f"  Exception: {exception}")
    
    def log_assembly_stats(self, fasta_path):
        """Log final assembly statistics."""
        try:
            from Bio import SeqIO
            
            contigs = list(SeqIO.parse(str(fasta_path), "fasta"))
            lengths = [len(r.seq) for r in contigs]
            lengths.sort(reverse=True)
            
            total_length = sum(lengths)
            n_contigs = len(lengths)
            
            # Calculate N50
            cumsum = 0
            n50 = 0
            for length in lengths:
                cumsum += length
                if cumsum >= total_length / 2:
                    n50 = length
                    break
            
            stats = {
                "n_contigs": n_contigs,
                "total_length": total_length,
                "n50": n50,
                "largest_contig": max(lengths) if lengths else 0,
                "smallest_contig": min(lengths) if lengths else 0
            }
            
            self.run_data["final_stats"]["assembly"] = stats
            self._append_log(f"ASSEMBLY STATS: {n_contigs} contigs, {total_length:,} bp, N50={n50:,}")
            
        except Exception as e:
            self.log_error(f"Failed to calculate assembly stats: {e}")
    
    def finalize_log(self):
        """Finalize the log with total runtime and summary."""
        end_time = time.time()
        total_runtime = end_time - self.start_time
        
        self.run_data["end_time"] = datetime.now().isoformat()
        self.run_data["total_runtime_seconds"] = total_runtime
        
        # Write summary
        summary = f"""

{'='*80}
RUN SUMMARY
{'='*80}
Total Runtime: {total_runtime/60:.1f} minutes
Commands Executed: {len(self.run_data['commands_executed'])}
Modules Completed: {len(self.run_data['module_results'])}
Errors Encountered: {len(self.run_data['errors'])}

Final Assembly Stats:
"""
        
        if "assembly" in self.run_data["final_stats"]:
            stats = self.run_data["final_stats"]["assembly"]
            summary += f"  Contigs: {stats['n_contigs']}\n"
            summary += f"  Total Length: {stats['total_length']:,} bp\n"
            summary += f"  N50: {stats['n50']:,} bp\n"
        
        summary += f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        summary += "="*80 + "\n"
        
        self._append_log(summary)
        
        # Write JSON data file for programmatic access
        json_file = self.outdir / "funcat_run_data.json"
        with open(json_file, 'w') as f:
            json.dump(self.run_data, f, indent=2)
    
    def _append_log(self, message):
        """Append message to log file."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(self.log_file, 'a') as f:
            f.write(f"[{timestamp}] {message}\n")


# Global logger instance
_logger = None

def init_logger(outdir):
    """Initialize global logger."""
    global _logger
    _logger = FuncatLogger(outdir)
    return _logger

def get_logger():
    """Get current logger instance."""
    return _logger

def log_user_selection(step, selection, value):
    """Log user selection (convenience function)."""
    if _logger:
        _logger.log_user_selection(step, selection, value)

def log_command(command, cwd=None):
    """Log command execution (convenience function)."""
    if _logger:
        _logger.log_command(command, cwd)

def log_module_start(module_name):
    """Log module start (convenience function)."""
    if _logger:
        _logger.log_module_start(module_name)

def log_module_end(module_name, result_path=None, stats=None):
    """Log module end (convenience function)."""
    if _logger:
        _logger.log_module_end(module_name, result_path, stats)

def log_error(error_msg, exception=None):
    """Log error (convenience function)."""
    if _logger:
        _logger.log_error(error_msg, exception)

def finalize_log():
    """Finalize log (convenience function)."""
    if _logger:
        _logger.finalize_log()
