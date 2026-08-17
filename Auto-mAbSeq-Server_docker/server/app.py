import os
import uuid
import shutil
import subprocess
import json
from datetime import datetime
from pathlib import Path
from threading import Thread
from flask import Flask, render_template, request, jsonify, send_file
from flask_cors import CORS
import yaml

app = Flask(__name__)
CORS(app)

# ===== PATHS (можно переопределить переменными окружения для Docker) =====
BASE_DIR = Path(__file__).parent
DATA_DIR = Path(os.environ.get('DATA_DIR', BASE_DIR / 'data'))
JOBS_DIR = DATA_DIR / 'jobs'
ARCHIVES_DIR = DATA_DIR / 'archives'
RAW_DATA_DIR = DATA_DIR / 'raw' / 'sanger'
LOGS_DIR = DATA_DIR / 'logs'
RESULTS_DIR = DATA_DIR / 'results'
INTERMEDIATE_DIR = DATA_DIR / 'intermediate'
CONFIG_FILE = BASE_DIR / 'config' / 'config.yaml'
PIPELINE_SCRIPT = BASE_DIR / 'run_full_pipeline.sh'

# ===== ENSURE DIRECTORIES EXIST =====
for dir_path in [JOBS_DIR, ARCHIVES_DIR, RAW_DATA_DIR, LOGS_DIR, RESULTS_DIR, INTERMEDIATE_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

TITLE = "="*80

def load_config(config_file=CONFIG_FILE):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def save_config(config, config_file=CONFIG_FILE):
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

def cleanup_pipeline_data():
    """Очистить общие директории (перед загрузкой новых файлов)"""
    try:
        if RAW_DATA_DIR.exists():
            for file in RAW_DATA_DIR.glob('*'):
                if file.is_file():
                    file.unlink()
                elif file.is_dir():
                    shutil.rmtree(file)
        
        if LOGS_DIR.exists():
            for file in LOGS_DIR.glob('*'):
                if file.is_file():
                    file.unlink()
                elif file.is_dir():
                    shutil.rmtree(file)
        
        if RESULTS_DIR.exists():
            for item in RESULTS_DIR.iterdir():
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
        
        if INTERMEDIATE_DIR.exists():
            for item in INTERMEDIATE_DIR.iterdir():
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
        
        app.logger.info("Pipeline data cleaned up successfully")
    except Exception as e:
        app.logger.error(f"Cleanup error: {str(e)}")

def get_job_dir(job_id):
    return JOBS_DIR / job_id

def get_job_status_file(job_id):
    return get_job_dir(job_id) / 'status.json'

def get_job_log_file(job_id):
    return get_job_dir(job_id) / 'pipeline.log'

def init_job(job_id):
    """Инициализировать job и создать status файл"""
    job_dir = get_job_dir(job_id)
    job_dir.mkdir(parents=True, exist_ok=True)
    
    # Создать директории для этого job'а
    (job_dir / 'data' / 'raw' / 'sanger').mkdir(parents=True, exist_ok=True)
    (job_dir / 'data' / 'intermediate').mkdir(parents=True, exist_ok=True)
    (job_dir / 'data' / 'results').mkdir(parents=True, exist_ok=True)
    (job_dir / 'data' / 'logs').mkdir(parents=True, exist_ok=True)
    
    status = {
        'job_id': job_id,
        'status': 'queued',
        'created_at': datetime.now().isoformat(),
        'started_at': None,
        'completed_at': None,
        'progress': 0,
        'error': None,
        'message': 'Waiting to process...'
    }
    
    with open(get_job_status_file(job_id), 'w') as f:
        json.dump(status, f, indent=2)
    
    return status

def update_job_status(job_id, status, progress=None, message=None, error=None):
    """Обновить статус job'а"""
    status_file = get_job_status_file(job_id)
    
    with open(status_file, 'r') as f:
        data = json.load(f)
    
    data['status'] = status
    
    if progress is not None:
        data['progress'] = progress
    
    if message is not None:
        data['message'] = message
    
    if error is not None:
        data['error'] = error
    
    if status == 'running' and not data['started_at']:
        data['started_at'] = datetime.now().isoformat()
    
    if status == 'completed' and not data['completed_at']:
        data['completed_at'] = datetime.now().isoformat()
    
    with open(status_file, 'w') as f:
        json.dump(data, f, indent=2)

def get_job_status(job_id):
    """Получить статус job'а"""
    status_file = get_job_status_file(job_id)
    
    if not status_file.exists():
        return None
    
    with open(status_file, 'r') as f:
        return json.load(f)

def create_archive(job_id):
    """Создать ZIP архив результатов"""
    job_dir = get_job_dir(job_id)
    archive_path = ARCHIVES_DIR / f'{job_id}.zip'
    
    try:
        shutil.make_archive(str(archive_path.with_suffix('')), 'zip', job_dir, 'data')
        return archive_path
    except Exception as e:
        raise Exception(f"Failed to create archive: {str(e)}")

def cleanup_job(job_id):
    """Удалить job директорию и архив"""
    job_dir = get_job_dir(job_id)
    archive_path = ARCHIVES_DIR / f'{job_id}.zip'
    
    try:
        if job_dir.exists():
            shutil.rmtree(job_dir)
        if archive_path.exists():
            archive_path.unlink()
    except Exception as e:
        app.logger.error(f"Cleanup error for job {job_id}: {str(e)}")

def run_pipeline(job_id, params=None):
    """
    Запустить пайплайн в ИЗОЛИРОВАННОЙ директории для этого job'а
    Использует Threading для неблокирования Flask
    """
    job_dir = get_job_dir(job_id)
    task_data_dir = job_dir / 'data'
    
    try:
        # Обновить конфиг если переданы параметры
        if params:
            config = load_config()
            
            if 'forward_line' in params:
                config['ab1_conversion']['forward_line'] = params['forward_line']
            
            if 'min_aa' in params:
                config['anarci']['min_aa'] = int(params['min_aa'])
            
            if 'max_aa' in params:
                config['anarci']['max_aa'] = int(params['max_aa'])
            
            if 'scheme' in params:
                config['anarci']['scheme'] = params['scheme']
            
            if 'tool' in params:
                config['alignment']['tool'] = params['tool']
            
            save_config(config)
        
        update_job_status(job_id, 'running', 0, 'Starting pipeline...')
        
        # ===== КЛЮЧЕВОЙ МОМЕНТ: Переменные окружения для изоляции =====
        log_file = get_job_log_file(job_id)
        
        # Создать окружение с переопределённые директориями
        env = os.environ.copy()
        env['DATA_DIR'] = str(task_data_dir)
        env['RAW_DATA_DIR'] = str(task_data_dir / 'raw' / 'sanger')
        env['INTERMEDIATE_DIR'] = str(task_data_dir / 'intermediate')
        env['RESULTS_DIR'] = str(task_data_dir / 'results')
        env['LOGS_DIR'] = str(task_data_dir / 'logs')
        env['JOB_ID'] = job_id
        
        # Исполнить пайплайн с переопределёнными путями
        with open(log_file, 'w') as log:
            result = subprocess.run(
                [str(PIPELINE_SCRIPT), 'config/config.yaml'],
                cwd=str(BASE_DIR),
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True
            )
        
        if result.returncode == 0:
            update_job_status(job_id, 'completed', 100, 'Pipeline completed successfully!')
            
            # Результаты уже в task_data_dir, ничего копировать не нужно
            app.logger.info(f"Job {job_id} completed successfully")
        else:
            with open(log_file, 'r') as log:
                error_msg = log.read()[-500:]
            
            update_job_status(job_id, 'failed', 0, 'Pipeline failed', error_msg)
            app.logger.error(f"Job {job_id} failed")
    
    except Exception as e:
        update_job_status(job_id, 'failed', 0, 'Pipeline error', str(e))
        app.logger.error(f"Job {job_id} error: {str(e)}")

def run_pipeline_threaded(job_id, params=None):
    """Запустить пайплайн в отдельном потоке"""
    thread = Thread(target=run_pipeline, args=(job_id, params), daemon=False)
    thread.start()

# ===== ROUTES =====

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/config', methods=['GET'])
def get_ui_config():
    """Получить UI конфиг с параметрами"""
    config = load_config()
    
    ui_config = {
        'parameters': {
            'forward_line': {
                'type': 'toggle',
                'label': 'Forward Strand',
                'default': config['ab1_conversion'].get('forward_line', False),
                'description': 'False = reverse complement'
            },
            'min_aa': {
                'type': 'slider',
                'min': 50,
                'max': 120,
                'default': config['anarci'].get('min_aa', 80),
                'label': 'Min AA length'
            },
            'max_aa': {
                'type': 'slider',
                'min': 100,
                'max': 200,
                'default': config['anarci'].get('max_aa', 150),
                'label': 'Max AA length'
            },
            'scheme': {
                'type': 'dropdown',
                'options': ['imgt', 'kabat', 'chothia'],
                'default': config['anarci'].get('scheme', 'imgt'),
                'label': 'Numbering scheme'
            },
            'tool': {
                'type': 'dropdown',
                'options': ['clustalo', 'mafft'],
                'default': config['alignment'].get('tool', 'clustalo'),
                'label': 'Alignment tool'
            }
        }
    }
    
    return jsonify(ui_config)

@app.route('/api/upload', methods=['POST'])
def upload_files():
    """Загрузить файлы и запустить пайплайн"""
    
    # Очистить общие директории перед новой загрузкой
    cleanup_pipeline_data()
    
    if 'files' not in request.files:
        return jsonify({'error': 'No files provided'}), 400
    
    files = request.files.getlist('files')
    
    if not files or len(files) == 0:
        return jsonify({'error': 'No files selected'}), 400
    
    if len(files) > 100:
        return jsonify({'error': 'Maximum 100 files allowed'}), 400
    
    # Инициализировать job
    job_id = str(uuid.uuid4())
    init_job(job_id)
    
    try:
        # Сохранить файлы в job-specific директорию
        job_raw_dir = get_job_dir(job_id) / 'data' / 'raw' / 'sanger'
        
        for file in files:
            if file.filename.endswith('.ab1'):
                file.save(job_raw_dir / file.filename)
        
        # Получить параметры из формы
        params = {
            'forward_line': request.form.get('forward_line', 'false').lower() == 'true',
            'min_aa': request.form.get('min_aa', '80'),
            'max_aa': request.form.get('max_aa', '150'),
            'scheme': request.form.get('scheme', 'imgt'),
            'tool': request.form.get('tool', 'clustalo')
        }
        
        # Запустить пайплайн в отдельном потоке
        run_pipeline_threaded(job_id, params)
        
        return jsonify({
            'job_id': job_id,
            'status': 'queued',
            'message': 'Files uploaded, pipeline starting...'
        }), 202
    
    except Exception as e:
        update_job_status(job_id, 'failed', 0, 'Upload error', str(e))
        return jsonify({'error': str(e)}), 500

@app.route('/api/jobs/<job_id>/status', methods=['GET'])
def job_status(job_id):
    """Получить статус job'а"""
    status = get_job_status(job_id)
    
    if not status:
        return jsonify({'error': 'Job not found'}), 404
    
    return jsonify(status)

@app.route('/api/jobs/<job_id>/download', methods=['GET'])
def download_results(job_id):
    """Скачать результаты"""
    status = get_job_status(job_id)
    
    if not status:
        return jsonify({'error': 'Job not found'}), 404
    
    if status['status'] != 'completed':
        return jsonify({'error': f"Job is {status['status']}, cannot download"}), 400
    
    try:
        # Создать архив
        archive_path = create_archive(job_id)
        
        # Отправить файл
        response = send_file(
            archive_path,
            as_attachment=True,
            download_name=f"mAbSeq_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
        )
        
        # Очистить job после скачивания
        cleanup_job(job_id)
        cleanup_pipeline_data()
        
        return response
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/jobs/<job_id>/delete', methods=['DELETE'])
def delete_job(job_id):
#    """Удалить job"""
    try:
        cleanup_job(job_id)
        cleanup_pipeline_data()
        return jsonify({'message': 'Job deleted'}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/health', methods=['GET'])
def healthcheck():
    """Health check"""
    return jsonify({'status': 'ok'}), 200

# ===== ERROR HANDLERS =====

@app.errorhandler(404)
def not_found_error(error):
    return jsonify({'error': 'Not found'}), 404

@app.errorhandler(500)
def internal_error(error):
    return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000, threaded=True)
